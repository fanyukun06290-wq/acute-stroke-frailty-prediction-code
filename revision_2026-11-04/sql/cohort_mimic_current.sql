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
WITH first_icu AS (
 SELECT i.*,row_number() OVER(PARTITION BY hadm_id ORDER BY intime,stay_id) AS icu_order
 FROM mimiciv_icu.icustays i
), dx_tokens AS (
    SELECT
        d.hadm_id,
        d.icd_version,
        replace(upper(trim(d.icd_code)), '.', '') AS code
    FROM mimiciv_hosp.diagnoses_icd AS d
), dx_flags AS (
    SELECT
        hadm_id,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code IN
                ('43301','43311','43321','43331','43381','43391',
                 '43401','43411','43491'))
            OR (icd_version = 10 AND code LIKE 'I63%'
                AND code NOT LIKE 'I636%')
        ), FALSE) AS has_ais,
        COALESCE(BOOL_OR(
            (icd_version = 10 AND code LIKE 'I636%')
            OR (icd_version = 9 AND code = '4376')
            OR (icd_version = 10 AND code LIKE 'I676%')
        ), FALSE) AS has_cvst,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code = '431')
            OR (icd_version = 10 AND code LIKE 'I61%')
        ), FALSE) AS has_ich,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code = '430')
            OR (icd_version = 10 AND code LIKE 'I60%')
        ), FALSE) AS has_sah,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code = '436')
            OR (icd_version = 10 AND code LIKE 'I64%')
        ), FALSE) AS has_unspecified_stroke,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '432%')
            OR (icd_version = 10 AND code LIKE 'I62%')
        ), FALSE) AS has_other_intracranial_hemorrhage,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '85%')
            OR (icd_version = 10 AND code LIKE 'S06%')
        ), FALSE) AS has_traumatic_intracranial_injury,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '435%')
            OR (icd_version = 10 AND code LIKE 'G45%')
        ), FALSE) AS has_tia,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '438%')
            OR (icd_version = 10 AND code LIKE 'I69%')
            OR (icd_version = 9 AND code = 'V1254')
            OR (icd_version = 10 AND code LIKE 'Z8673%')
        ), FALSE) AS has_history_or_sequela
FROM dx_tokens GROUP BY hadm_id
), cohort_base AS (
 SELECT 'MIMIC-IV'::text AS data_source,i.subject_id,i.hadm_id,i.stay_id,
 'MIMIC:'||i.subject_id::text AS patient_group_id,
 i.intime AS icu_intime,i.outtime AS icu_outtime,
 a.admittime,a.dischtime,a.deathtime,a.hospital_expire_flag,
 upper(trim(p.gender)) AS sex,
 (p.anchor_age+extract(year FROM a.admittime)::integer-p.anchor_year)::numeric AS age_years,
 COALESCE(f.has_ais,false) AS has_ais,
COALESCE(f.has_ich,false) AS has_ich,
COALESCE(f.has_sah,false) AS has_sah,
COALESCE(f.has_cvst,false) AS has_cvst,
COALESCE(f.has_unspecified_stroke,false) AS has_unspecified_stroke,
COALESCE(f.has_other_intracranial_hemorrhage,false) AS has_other_intracranial_hemorrhage,
COALESCE(f.has_traumatic_intracranial_injury,false) AS has_traumatic_intracranial_injury
 FROM first_icu i JOIN mimiciv_hosp.admissions a USING(hadm_id)
 JOIN mimiciv_hosp.patients p ON p.subject_id=i.subject_id
 LEFT JOIN dx_flags f ON f.hadm_id=i.hadm_id WHERE i.icu_order=1
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
 (icu_intime IS NULL OR icu_outtime IS NULL OR admittime IS NULL OR dischtime IS NULL
  OR admittime>icu_intime OR icu_outtime<=icu_intime OR dischtime<=icu_intime
  OR deathtime<icu_intime) IS TRUE AS invalid_time_flag,
 CASE WHEN hospital_expire_flag NOT IN (0,1) OR hospital_expire_flag IS NULL
  OR dischtime IS NULL OR dischtime<=icu_intime
  OR (hospital_expire_flag=1 AND (deathtime IS NULL OR deathtime<icu_intime OR deathtime>dischtime))
  OR (hospital_expire_flag=0 AND deathtime IS NOT NULL)
 THEN 1 ELSE 0 END AS outcome_missing_flag,
 LEAST(icu_outtime,dischtime,deathtime) AS effective_icu_end_at,
 extract(epoch FROM (deathtime-icu_intime))/60.0 AS death_offset_min,
 extract(epoch FROM (dischtime-icu_intime))/60.0 AS hospital_discharge_offset_min,
 (icu_outtime>dischtime)::integer AS icu_end_after_hospital_flag
 FROM labeled WHERE n_primary_acute_groups>0 AND NOT has_traumatic_intracranial_injury
 AND age_years>=18
), primary_cohort AS (
 SELECT *,extract(epoch FROM (effective_icu_end_at-icu_intime))/60.0 AS effective_icu_end_min,
 CASE WHEN NOT invalid_time_flag AND effective_icu_end_at>icu_intime+interval '6 hours'
 THEN 1 ELSE 0 END AS landmark_6h_eligible,
 CASE WHEN NOT invalid_time_flag AND effective_icu_end_at>icu_intime+interval '24 hours'
 THEN 1 ELSE 0 END AS landmark_24h_eligible,
 CASE WHEN outcome_missing_flag=1 THEN NULL
 ELSE (hospital_expire_flag=1 AND death_offset_min>360 AND death_offset_min<43200)::integer END
 AS death_6h_to_30d_inhospital,
 CASE WHEN outcome_missing_flag=1 THEN NULL
 ELSE (hospital_expire_flag=1 AND death_offset_min>1440 AND death_offset_min<43200)::integer END
 AS death_24h_to_30d_inhospital
 FROM timed
) SELECT * FROM primary_cohort ORDER BY stay_id;
