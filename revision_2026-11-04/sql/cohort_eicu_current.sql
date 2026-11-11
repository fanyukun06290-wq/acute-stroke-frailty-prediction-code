/* Reviewer revision route A, confirmed by author 2026-11-04.
READ ONLY. Do not use phenotype or metadata columns as model predictors.
Index stay: first ICU/unit in each hospital admission, before stroke selection.
Landmark: alive and still in BOTH ICU/unit and hospital strictly after hour 6/24.
Predictors: measurements in [0, landmark), with known recording time also < landmark.
Missing recording time: unavailable for timestamped lab/chartevent sources.
eICU monitor vitals: observation time only; separate entry time is not supplied.
Endpoint: recorded in-hospital death strictly after landmark and before ICU day 30.
No comorbidity variables. Stroke phenotype is retrospective cohort metadata only.
Unknown/contradictory outcomes are flagged, never silently converted to survivors.
Output includes source rows with eligibility flags for auditable cohort flow.
*/
WITH first_unit AS (
 SELECT p.*,row_number() OVER(PARTITION BY patienthealthsystemstayid
 ORDER BY hospitaladmitoffset DESC NULLS LAST,unitvisitnumber NULLS LAST,patientunitstayid) AS unit_order,
 count(*) FILTER(WHERE hospitaladmitoffset IS NULL) OVER(PARTITION BY patienthealthsystemstayid)
 AS missing_order_offsets,
 count(*) OVER(PARTITION BY patienthealthsystemstayid,hospitaladmitoffset) AS same_offset_units
 FROM eicu_crd.patient p
), eligible_units AS (
    SELECT patientunitstayid
    FROM first_unit
    WHERE unit_order = 1
), tokenized AS (
    SELECT
        d.patientunitstayid,
        d.diagnosisoffset,
        d.diagnosisstring,
        upper(trim(tok.raw_token)) AS raw_token,
        replace(upper(trim(tok.raw_token)), '.', '') AS code,
        lower(coalesce(d.diagnosisstring, '')) AS dx_text
    FROM eicu_crd.diagnosis AS d
    INNER JOIN eligible_units AS eu
        ON eu.patientunitstayid = d.patientunitstayid
    CROSS JOIN LATERAL regexp_split_to_table(
        coalesce(d.icd9code, ''),
        '[[:space:]]*,[[:space:]]*'
    ) AS tok(raw_token)
    WHERE nullif(trim(tok.raw_token), '') IS NOT NULL
), classified AS (
    SELECT
        *,
        (code IN
            ('43301','43311','43321','43331','43381','43391',
             '43401','43411','43491')) AS is_ais_icd9,
        (code = '431'
         AND dx_text NOT LIKE '%trauma - cns%'
         AND dx_text NOT LIKE '%traumatic%'
         AND dx_text NOT LIKE '%intracranial injury%') AS is_ich_icd9,
        (code = '430'
         AND dx_text NOT LIKE '%trauma - cns%'
         AND dx_text NOT LIKE '%traumatic%'
         AND dx_text NOT LIKE '%intracranial injury%') AS is_sah_icd9,
        (code = '436') AS is_unspecified_icd9,
        (code LIKE '432%'
         AND dx_text NOT LIKE '%trauma - cns%'
         AND dx_text NOT LIKE '%traumatic%'
         AND dx_text NOT LIKE '%intracranial injury%') AS is_other_ich_icd9,
        (code LIKE 'I63%' AND code NOT LIKE 'I636%')
            AS is_ais_icd10_concordant,
        (code LIKE 'I61%'
         AND dx_text NOT LIKE '%trauma - cns%'
         AND dx_text NOT LIKE '%traumatic%'
         AND dx_text NOT LIKE '%intracranial injury%')
            AS is_ich_icd10_concordant,
        (code LIKE 'I60%'
         AND dx_text NOT LIKE '%trauma - cns%'
         AND dx_text NOT LIKE '%traumatic%'
         AND dx_text NOT LIKE '%intracranial injury%')
            AS is_sah_icd10_concordant,
        (code LIKE 'I64%') AS is_unspecified_icd10_concordant,
        (code LIKE 'I62%'
         AND dx_text NOT LIKE '%trauma - cns%'
         AND dx_text NOT LIKE '%traumatic%'
         AND dx_text NOT LIKE '%intracranial injury%')
            AS is_other_ich_icd10_concordant,
        (code = '4376'
         OR code LIKE 'I636%'
         OR code LIKE 'I676%'
         OR dx_text LIKE '%cortical venous thrombosis%'
         OR dx_text LIKE '%cerebral venous thrombosis%'
         OR dx_text LIKE '%intracranial venous thrombosis%'
         OR dx_text LIKE '%thrombosis of intracranial venous%') AS is_cvst,
        (code LIKE '85%'
         OR code LIKE 'S06%'
         OR dx_text LIKE '%trauma - cns%'
         OR dx_text LIKE '%traumatic%'
         OR dx_text LIKE '%intracranial injury%') AS is_trauma
    FROM tokenized
), flags AS (
    SELECT
        patientunitstayid,
        COALESCE(BOOL_OR(is_ais_icd9), FALSE) AS has_ais,
        COALESCE(BOOL_OR(is_ich_icd9), FALSE) AS has_ich,
        COALESCE(BOOL_OR(is_sah_icd9), FALSE) AS has_sah,
        COALESCE(BOOL_OR(is_cvst), FALSE) AS has_cvst,
        COALESCE(BOOL_OR(is_unspecified_icd9), FALSE)
            AS has_unspecified_stroke,
        COALESCE(BOOL_OR(is_other_ich_icd9), FALSE)
            AS has_other_intracranial_hemorrhage,
        COALESCE(BOOL_OR(is_trauma), FALSE)
            AS has_traumatic_intracranial_injury,
        COALESCE(BOOL_OR(is_ais_icd10_concordant), FALSE)
            AS has_ais_icd10_concordant,
        COALESCE(BOOL_OR(is_ich_icd10_concordant), FALSE)
            AS has_ich_icd10_concordant,
        COALESCE(BOOL_OR(is_sah_icd10_concordant), FALSE)
            AS has_sah_icd10_concordant,
        MIN(diagnosisoffset) FILTER (WHERE
            is_ais_icd9 OR is_ich_icd9 OR is_sah_icd9
            OR is_unspecified_icd9 OR is_other_ich_icd9 OR is_cvst
            OR is_ais_icd10_concordant OR is_ich_icd10_concordant
            OR is_sah_icd10_concordant) AS first_qualifying_code_offset
    FROM classified
    GROUP BY patientunitstayid
), cohort_base AS (
 SELECT 'eICU-CRD'::text AS data_source,p.patientunitstayid AS stay_id,p.patientunitstayid,p.patienthealthsystemstayid,
 p.uniquepid,'eICU:'||p.uniquepid AS patient_group_id,p.hospitalid,p.unitvisitnumber,p.unittype,
 p.hospitaladmitoffset,p.hospitaldischargeoffset,p.unitdischargeoffset,
 lower(trim(p.hospitaldischargestatus)) AS hospital_status,
 lower(trim(p.unitdischargestatus)) AS unit_status,
 p.missing_order_offsets,p.same_offset_units,
 CASE WHEN trim(p.age)~'^[0-9]+$' THEN trim(p.age)::numeric
 WHEN trim(p.age) LIKE '>%' THEN 90::numeric ELSE NULL::numeric END AS age_years,
 CASE WHEN upper(trim(p.gender)) IN ('MALE','M') THEN 'M'
 WHEN upper(trim(p.gender)) IN ('FEMALE','F') THEN 'F' ELSE NULL END AS sex,
 COALESCE(f.has_ais,false) AS has_ais,
COALESCE(f.has_ich,false) AS has_ich,
COALESCE(f.has_sah,false) AS has_sah,
COALESCE(f.has_cvst,false) AS has_cvst,
COALESCE(f.has_unspecified_stroke,false) AS has_unspecified_stroke,
COALESCE(f.has_other_intracranial_hemorrhage,false) AS has_other_intracranial_hemorrhage,
COALESCE(f.has_traumatic_intracranial_injury,false) AS has_traumatic_intracranial_injury,f.first_qualifying_code_offset
 FROM first_unit p LEFT JOIN flags f USING(patientunitstayid) WHERE p.unit_order=1
), derived AS (
    SELECT
        *,
        CAST(has_ais AS integer)
        + CAST(has_ich AS integer)
        + CAST(has_sah AS integer) AS n_primary_acute_groups
    FROM cohort_base
), labeled AS (
    SELECT
        *,
        CASE
            WHEN has_traumatic_intracranial_injury
                 AND (n_primary_acute_groups > 0 OR has_cvst
                      OR has_unspecified_stroke
                      OR has_other_intracranial_hemorrhage)
                THEN 'EXCLUDE_TRAUMATIC'
            WHEN n_primary_acute_groups > 1 THEN 'MIXED_PRIMARY_ACUTE'
            WHEN n_primary_acute_groups = 1 AND has_cvst
                THEN 'MIXED_WITH_CVST'
            WHEN n_primary_acute_groups = 1
                 AND (has_unspecified_stroke
                      OR has_other_intracranial_hemorrhage)
                THEN 'MIXED_WITH_AMBIGUOUS_SECONDARY'
            WHEN has_cvst THEN 'CVST'
            WHEN has_ais THEN 'AIS_STRICT'
            WHEN has_ich THEN 'ICH_STRICT'
            WHEN has_sah THEN 'SAH'
            WHEN has_unspecified_stroke THEN 'UNSPECIFIED_STROKE'
            WHEN has_other_intracranial_hemorrhage
                THEN 'OTHER_INTRACRANIAL_HEMORRHAGE'
            ELSE 'NO_STROKE_GROUP'
        END AS audit_phenotype
    FROM derived
), timed AS (
 SELECT *,
 (hospitaladmitoffset IS NULL OR hospitaladmitoffset>0 OR unitdischargeoffset IS NULL
 OR unitdischargeoffset<=0 OR hospitaldischargeoffset IS NULL OR hospitaldischargeoffset<=0
 OR missing_order_offsets>0) AS invalid_time_flag,
 LEAST(unitdischargeoffset,hospitaldischargeoffset)::numeric AS effective_icu_end_min,
 hospitaldischargeoffset::numeric AS hospital_discharge_offset_min,
 (unitdischargeoffset>hospitaldischargeoffset)::integer AS icu_end_after_hospital_flag,
 CASE WHEN hospital_status='expired' THEN
   CASE WHEN unit_status='expired' THEN LEAST(unitdischargeoffset,hospitaldischargeoffset)
   ELSE hospitaldischargeoffset END ELSE NULL END::numeric AS death_offset_min,
 CASE WHEN hospital_status='expired' AND unit_status='expired' THEN 'unit_death_discharge'
 WHEN hospital_status='expired' THEN 'hospital_death_discharge' ELSE NULL END AS death_time_source,
 CASE WHEN hospital_status IS NULL OR hospital_status NOT IN ('alive','expired')
 OR hospitaldischargeoffset IS NULL OR hospitaldischargeoffset<=0
 OR (unit_status='expired' AND hospital_status<>'expired') THEN 1 ELSE 0 END AS outcome_missing_flag
 FROM labeled WHERE n_primary_acute_groups>0 AND NOT has_traumatic_intracranial_injury
 AND age_years>=18
), primary_cohort AS (
 SELECT *,
 CASE WHEN NOT invalid_time_flag AND effective_icu_end_min>360 THEN 1 ELSE 0 END AS landmark_6h_eligible,
 CASE WHEN NOT invalid_time_flag AND effective_icu_end_min>1440 THEN 1 ELSE 0 END AS landmark_24h_eligible,
 CASE WHEN outcome_missing_flag=1 THEN NULL
 ELSE (hospital_status='expired' AND death_offset_min>360 AND death_offset_min<43200)::integer END
 AS death_6h_to_30d_inhospital,
 CASE WHEN outcome_missing_flag=1 THEN NULL
 ELSE (hospital_status='expired' AND death_offset_min>1440 AND death_offset_min<43200)::integer END
 AS death_24h_to_30d_inhospital FROM timed
) SELECT * FROM primary_cohort ORDER BY stay_id;
