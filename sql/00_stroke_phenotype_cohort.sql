/*
FINAL ANALYSIS COHORT: ICU acute stroke prognosis
Databases: MIMIC-IV v3.1 and eICU-CRD v2.0

Purpose:
  Build one row per first ICU stay within each hospital admission/stay.
  The result contains the primary analysis cohort plus phenotype flags,
  subtype eligibility, diagnosis timing, and 30-day in-hospital mortality.

Primary cohort rule:
  - Include strict AIS, ICH, and SAH.
  - Retain mixed/ambiguous acute stroke records in the overall cohort,
    but mark them as ineligible for mutually exclusive subtype analyses.
  - Exclude traumatic intracranial injury/hemorrhage from the primary cohort.
  - Exclude CVST, unspecified stroke, and other intracranial hemorrhage from
    the primary cohort when they occur without a strict AIS/ICH/SAH code.
  - Those groups remain identifiable for prespecified sensitivity analyses.

Primary outcome:
  30-day in-hospital mortality from ICU admission.
  This is not 30-day post-discharge all-cause mortality.

Important execution rule:
  Run the MIMIC block against the MIMIC-IV v3.1 schema and the eICU block
  against the eICU-CRD v2.0 schema. Do not execute both blocks on one connection.

This file is READ-ONLY. It contains no CREATE, INSERT, UPDATE, DELETE,
TRUNCATE, ALTER, or DROP statement.
*/


/* ================================================================
   BLOCK M. MIMIC-IV v3.1 final analysis cohort
   ================================================================ */

WITH first_icu AS (
    SELECT
        i.subject_id,
        i.hadm_id,
        i.stay_id,
        i.intime AS icu_intime,
        ROW_NUMBER() OVER (
            PARTITION BY i.hadm_id
            ORDER BY i.intime ASC, i.stay_id
        ) AS icu_order
    FROM mimiciv_icu.icustays AS i
), dx_tokens AS (
    SELECT
        d.hadm_id,
        d.icd_version,
        replace(upper(trim(d.icd_code)), '.', '') AS code
    FROM mimiciv_hosp.diagnoses_icd AS d
), dx_flags AS (
    SELECT
        hadm_id,

        /* Strict arterial AIS; ICD-9 436 is deliberately excluded. */
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code IN
                ('43301','43311','43321','43331','43381','43391',
                 '43401','43411','43491'))
            OR
            (icd_version = 10 AND code LIKE 'I63%'
                 AND code NOT LIKE 'I636%')
        ), FALSE) AS has_ais,

        /* Cerebral venous thrombosis is not arterial AIS. */
        COALESCE(BOOL_OR(
            icd_version = 10 AND code LIKE 'I636%'
        ), FALSE) AS has_cvst,

        /* Strict nontraumatic intracerebral hemorrhage. */
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code = '431')
            OR (icd_version = 10 AND code LIKE 'I61%')
        ), FALSE) AS has_ich,

        /* Strict nontraumatic subarachnoid hemorrhage. */
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code = '430')
            OR (icd_version = 10 AND code LIKE 'I60%')
        ), FALSE) AS has_sah,

        /* Ambiguous/sensitivity-only cerebrovascular groups. */
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code = '436')
            OR (icd_version = 10 AND code LIKE 'I64%')
        ), FALSE) AS has_unspecified_stroke,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '432%')
            OR (icd_version = 10 AND code LIKE 'I62%')
        ), FALSE) AS has_other_intracranial_hemorrhage,

        /* Non-acute/history codes are retained as flags, not exclusions. */
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '435%')
            OR (icd_version = 10 AND code LIKE 'G45%')
        ), FALSE) AS has_tia,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '438%')
            OR (icd_version = 10 AND code LIKE 'I69%')
            OR (icd_version = 9 AND code = 'V1254')
            OR (icd_version = 10 AND code LIKE 'Z8673%')
        ), FALSE) AS has_history_or_sequela,

        /* Conservative admission-level traumatic intracranial injury flag. */
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '85%')
            OR (icd_version = 10 AND code LIKE 'S06%')
        ), FALSE) AS has_traumatic_intracranial_injury
    FROM dx_tokens
    GROUP BY hadm_id
), cohort_base AS (
    SELECT
        i.subject_id,
        i.hadm_id,
        i.stay_id,
        i.icu_intime,
        a.admittime,
        a.dischtime,
        a.deathtime,
        a.hospital_expire_flag,
        p.anchor_year_group,
        COALESCE(f.has_ais, FALSE) AS has_ais,
        COALESCE(f.has_cvst, FALSE) AS has_cvst,
        COALESCE(f.has_ich, FALSE) AS has_ich,
        COALESCE(f.has_sah, FALSE) AS has_sah,
        COALESCE(f.has_unspecified_stroke, FALSE)
            AS has_unspecified_stroke,
        COALESCE(f.has_other_intracranial_hemorrhage, FALSE)
            AS has_other_intracranial_hemorrhage,
        COALESCE(f.has_tia, FALSE) AS has_tia,
        COALESCE(f.has_history_or_sequela, FALSE)
            AS has_history_or_sequela,
        COALESCE(f.has_traumatic_intracranial_injury, FALSE)
            AS has_traumatic_intracranial_injury
    FROM first_icu AS i
    INNER JOIN mimiciv_hosp.admissions AS a
        ON a.hadm_id = i.hadm_id
    INNER JOIN mimiciv_hosp.patients AS p
        ON p.subject_id = i.subject_id
    LEFT JOIN dx_flags AS f
        ON f.hadm_id = i.hadm_id
    WHERE i.icu_order = 1
), derived AS (
    SELECT
        *,
        (
            CAST(has_ais AS INTEGER)
            + CAST(has_ich AS INTEGER)
            + CAST(has_sah AS INTEGER)
        ) AS n_primary_acute_groups
    FROM cohort_base
), labeled AS (
    SELECT
        *,
        CASE
            WHEN has_traumatic_intracranial_injury
                 AND (
                     n_primary_acute_groups > 0
                     OR has_cvst
                     OR has_unspecified_stroke
                     OR has_other_intracranial_hemorrhage
                 )
                THEN 'EXCLUDE_TRAUMATIC'
            WHEN n_primary_acute_groups > 1
                THEN 'MIXED_PRIMARY_ACUTE'
            WHEN n_primary_acute_groups = 1 AND has_cvst
                THEN 'MIXED_WITH_CVST'
            WHEN n_primary_acute_groups = 1
                 AND (
                     has_unspecified_stroke
                     OR has_other_intracranial_hemorrhage
                 )
                THEN 'MIXED_WITH_AMBIGUOUS_SECONDARY'
            WHEN has_cvst
                THEN 'CVST'
            WHEN has_ais
                THEN 'AIS_STRICT'
            WHEN has_ich
                THEN 'ICH_STRICT'
            WHEN has_sah
                THEN 'SAH'
            WHEN has_unspecified_stroke
                THEN 'UNSPECIFIED_STROKE'
            WHEN has_other_intracranial_hemorrhage
                THEN 'OTHER_INTRACRANIAL_HEMORRHAGE'
            ELSE 'NO_STROKE_GROUP'
        END AS audit_phenotype
    FROM derived
), eligible AS (
    SELECT
        *,
        (
            n_primary_acute_groups > 0
            AND NOT has_traumatic_intracranial_injury
        ) AS primary_cohort_flag,
        (
            n_primary_acute_groups = 1
            AND NOT has_traumatic_intracranial_injury
            AND NOT has_cvst
            AND NOT has_unspecified_stroke
            AND NOT has_other_intracranial_hemorrhage
        ) AS subtype_eligible_flag
    FROM labeled
), final_cohort AS (
    SELECT
        'MIMIC-IV'::text AS data_source,
        subject_id,
        hadm_id,
        stay_id,
        anchor_year_group,
        icu_intime,
        admittime,
        dischtime,
        CASE
            WHEN dischtime >= icu_intime + INTERVAL '24 hours'
                THEN 1 ELSE 0
        END AS landmark_24h_eligible,
        audit_phenotype,
        primary_cohort_flag,
        subtype_eligible_flag,
        has_ais,
        has_ich,
        has_sah,
        has_cvst,
        has_unspecified_stroke,
        has_other_intracranial_hemorrhage,
        has_tia,
        has_history_or_sequela,
        has_traumatic_intracranial_injury,
        n_primary_acute_groups,
        hospital_expire_flag AS hospital_death_flag,
        CASE
            WHEN hospital_expire_flag = 1
                 AND deathtime IS NOT NULL
                 AND deathtime >= icu_intime
                 AND deathtime < icu_intime + INTERVAL '30 days'
                THEN 1 ELSE 0
        END AS death_30d_inhospital,
        CASE
            WHEN hospital_expire_flag = 1
                 AND deathtime IS NOT NULL
                 AND deathtime >= icu_intime + INTERVAL '24 hours'
                 AND deathtime < icu_intime + INTERVAL '30 days'
                THEN 1 ELSE 0
        END AS death_24h_to_30d_inhospital,
        CASE
            WHEN hospital_expire_flag IS NULL
                 OR (hospital_expire_flag = 1 AND deathtime IS NULL)
                THEN 1 ELSE 0
        END AS outcome_missing_flag
    FROM eligible
)
SELECT *
FROM final_cohort
WHERE primary_cohort_flag
ORDER BY hadm_id, stay_id;


/* ================================================================
   BLOCK E. eICU-CRD v2.0 final analysis cohort
   ================================================================ */

WITH first_unit AS (
    SELECT
        p.patientunitstayid,
        p.patienthealthsystemstayid,
        p.hospitalid,
        p.hospitaladmitoffset,
        p.hospitaldischargeoffset,
        p.hospitaldischargestatus,
        ROW_NUMBER() OVER (
            PARTITION BY p.patienthealthsystemstayid
            ORDER BY p.hospitaladmitoffset ASC NULLS LAST,
                     p.patientunitstayid
        ) AS unit_order
    FROM eicu_crd.patient AS p
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

        /* eICU primary definition uses the ICD-9 anchor token. */
        (
            code IN
                ('43301','43311','43321','43331','43381','43391',
                 '43401','43411','43491')
        ) AS is_ais_icd9,
        (
            code = '431'
            AND dx_text NOT LIKE '%trauma - cns%'
            AND dx_text NOT LIKE '%traumatic%'
            AND dx_text NOT LIKE '%intracranial injury%'
        ) AS is_ich_icd9,
        (
            code = '430'
            AND dx_text NOT LIKE '%trauma - cns%'
            AND dx_text NOT LIKE '%traumatic%'
            AND dx_text NOT LIKE '%intracranial injury%'
        ) AS is_sah_icd9,
        (code = '436') AS is_unspecified_icd9,
        (
            code LIKE '432%'
            AND dx_text NOT LIKE '%trauma - cns%'
            AND dx_text NOT LIKE '%traumatic%'
            AND dx_text NOT LIKE '%intracranial injury%'
        ) AS is_other_ich_icd9,

        /* ICD-10-looking tokens are concordance/sensitivity flags only. */
        (
            code LIKE 'I63%' AND code NOT LIKE 'I636%'
        ) AS is_ais_icd10_concordant,
        (
            code LIKE 'I61%'
            AND dx_text NOT LIKE '%trauma - cns%'
            AND dx_text NOT LIKE '%traumatic%'
            AND dx_text NOT LIKE '%intracranial injury%'
        ) AS is_ich_icd10_concordant,
        (
            code LIKE 'I60%'
            AND dx_text NOT LIKE '%trauma - cns%'
            AND dx_text NOT LIKE '%traumatic%'
            AND dx_text NOT LIKE '%intracranial injury%'
        ) AS is_sah_icd10_concordant,
        (code LIKE 'I64%') AS is_unspecified_icd10_concordant,
        (
            code LIKE 'I62%'
            AND dx_text NOT LIKE '%trauma - cns%'
            AND dx_text NOT LIKE '%traumatic%'
            AND dx_text NOT LIKE '%intracranial injury%'
        ) AS is_other_ich_icd10_concordant,

        /* CVST is not included in arterial AIS. */
        (
            code = '4376'
            OR code LIKE 'I636%'
            OR code LIKE 'I676%'
            OR dx_text LIKE '%cortical venous thrombosis%'
            OR dx_text LIKE '%cerebral venous thrombosis%'
            OR dx_text LIKE '%intracranial venous thrombosis%'
            OR dx_text LIKE '%thrombosis of intracranial venous%'
        ) AS is_cvst_related,

        /* Any explicit traumatic CNS signal excludes the primary cohort. */
        (
            code LIKE '85%'
            OR code LIKE 'S06%'
            OR dx_text LIKE '%trauma - cns%'
            OR dx_text LIKE '%traumatic%'
            OR dx_text LIKE '%intracranial injury%'
        ) AS is_traumatic_intracranial_signal
    FROM tokenized
), flags AS (
    SELECT
        patientunitstayid,
        COALESCE(BOOL_OR(is_ais_icd9), FALSE) AS has_ais,
        COALESCE(BOOL_OR(is_ich_icd9), FALSE) AS has_ich,
        COALESCE(BOOL_OR(is_sah_icd9), FALSE) AS has_sah,
        COALESCE(BOOL_OR(is_cvst_related), FALSE) AS has_cvst,
        COALESCE(BOOL_OR(is_unspecified_icd9), FALSE)
            AS has_unspecified_stroke,
        COALESCE(BOOL_OR(is_other_ich_icd9), FALSE)
            AS has_other_intracranial_hemorrhage,
        COALESCE(BOOL_OR(is_traumatic_intracranial_signal), FALSE)
            AS has_traumatic_intracranial_injury,

        /* Concordance flags for local mixed-token eICU data. */
        COALESCE(BOOL_OR(is_ais_icd10_concordant), FALSE)
            AS has_ais_icd10_concordant,
        COALESCE(BOOL_OR(is_ich_icd10_concordant), FALSE)
            AS has_ich_icd10_concordant,
        COALESCE(BOOL_OR(is_sah_icd10_concordant), FALSE)
            AS has_sah_icd10_concordant,
        COALESCE(BOOL_OR(is_unspecified_icd10_concordant), FALSE)
            AS has_unspecified_icd10_concordant,
        COALESCE(BOOL_OR(is_other_ich_icd10_concordant), FALSE)
            AS has_other_ich_icd10_concordant,

        MIN(diagnosisoffset) FILTER (WHERE
            is_ais_icd9 OR is_ich_icd9 OR is_sah_icd9
            OR is_unspecified_icd9 OR is_other_ich_icd9
            OR is_cvst_related
            OR is_ais_icd10_concordant
            OR is_ich_icd10_concordant
            OR is_sah_icd10_concordant
            OR is_unspecified_icd10_concordant
            OR is_other_ich_icd10_concordant
        ) AS first_qualifying_code_offset,
        MIN(diagnosisoffset) FILTER (WHERE
            is_ais_icd9 OR is_ich_icd9 OR is_sah_icd9
            OR is_unspecified_icd9 OR is_other_ich_icd9
        ) AS first_icd9_stroke_code_offset,
        MAX(CASE WHEN (
            is_ais_icd9 OR is_ich_icd9 OR is_sah_icd9
            OR is_unspecified_icd9 OR is_other_ich_icd9
            OR is_cvst_related
            OR is_ais_icd10_concordant
            OR is_ich_icd10_concordant
            OR is_sah_icd10_concordant
            OR is_unspecified_icd10_concordant
            OR is_other_ich_icd10_concordant
        ) AND diagnosisoffset >= 0 AND diagnosisoffset < 1440
            THEN 1 ELSE 0 END) AS qualifying_code_in_first_24h,
        MAX(CASE WHEN (
            is_ais_icd9 OR is_ich_icd9 OR is_sah_icd9
            OR is_unspecified_icd9 OR is_other_ich_icd9
        ) AND diagnosisoffset >= 0 AND diagnosisoffset < 1440
            THEN 1 ELSE 0 END) AS primary_icd9_code_in_first_24h
    FROM classified
    GROUP BY patientunitstayid
), cohort_base AS (
    SELECT
        p.patientunitstayid,
        p.patienthealthsystemstayid,
        p.hospitalid,
        p.hospitaladmitoffset,
        p.hospitaldischargeoffset,
        p.hospitaldischargestatus,
        COALESCE(f.has_ais, FALSE) AS has_ais,
        COALESCE(f.has_ich, FALSE) AS has_ich,
        COALESCE(f.has_sah, FALSE) AS has_sah,
        COALESCE(f.has_cvst, FALSE) AS has_cvst,
        COALESCE(f.has_unspecified_stroke, FALSE)
            AS has_unspecified_stroke,
        COALESCE(f.has_other_intracranial_hemorrhage, FALSE)
            AS has_other_intracranial_hemorrhage,
        COALESCE(f.has_traumatic_intracranial_injury, FALSE)
            AS has_traumatic_intracranial_injury,
        COALESCE(f.has_ais_icd10_concordant, FALSE)
            AS has_ais_icd10_concordant,
        COALESCE(f.has_ich_icd10_concordant, FALSE)
            AS has_ich_icd10_concordant,
        COALESCE(f.has_sah_icd10_concordant, FALSE)
            AS has_sah_icd10_concordant,
        COALESCE(f.has_unspecified_icd10_concordant, FALSE)
            AS has_unspecified_icd10_concordant,
        COALESCE(f.has_other_ich_icd10_concordant, FALSE)
            AS has_other_ich_icd10_concordant,
        f.first_qualifying_code_offset,
        f.first_icd9_stroke_code_offset,
        COALESCE(f.qualifying_code_in_first_24h, 0)
            AS qualifying_code_in_first_24h,
        COALESCE(f.primary_icd9_code_in_first_24h, 0)
            AS primary_icd9_code_in_first_24h
    FROM first_unit AS p
    LEFT JOIN flags AS f
        ON f.patientunitstayid = p.patientunitstayid
    WHERE p.unit_order = 1
), derived AS (
    SELECT
        *,
        (
            CAST(has_ais AS INTEGER)
            + CAST(has_ich AS INTEGER)
            + CAST(has_sah AS INTEGER)
        ) AS n_primary_acute_groups
    FROM cohort_base
), labeled AS (
    SELECT
        *,
        CASE
            WHEN has_traumatic_intracranial_injury
                 AND (
                     n_primary_acute_groups > 0
                     OR has_cvst
                     OR has_unspecified_stroke
                     OR has_other_intracranial_hemorrhage
                 )
                THEN 'EXCLUDE_TRAUMATIC'
            WHEN n_primary_acute_groups > 1
                THEN 'MIXED_PRIMARY_ACUTE'
            WHEN n_primary_acute_groups = 1 AND has_cvst
                THEN 'MIXED_WITH_CVST'
            WHEN n_primary_acute_groups = 1
                 AND (
                     has_unspecified_stroke
                     OR has_other_intracranial_hemorrhage
                 )
                THEN 'MIXED_WITH_AMBIGUOUS_SECONDARY'
            WHEN has_cvst
                THEN 'CVST'
            WHEN has_ais
                THEN 'AIS_STRICT'
            WHEN has_ich
                THEN 'ICH_STRICT'
            WHEN has_sah
                THEN 'SAH'
            WHEN has_unspecified_stroke
                THEN 'UNSPECIFIED_STROKE'
            WHEN has_other_intracranial_hemorrhage
                THEN 'OTHER_INTRACRANIAL_HEMORRHAGE'
            ELSE 'NO_STROKE_GROUP'
        END AS audit_phenotype
    FROM derived
), eligible AS (
    SELECT
        *,
        (
            n_primary_acute_groups > 0
            AND NOT has_traumatic_intracranial_injury
        ) AS primary_cohort_flag,
        (
            n_primary_acute_groups = 1
            AND NOT has_traumatic_intracranial_injury
            AND NOT has_cvst
            AND NOT has_unspecified_stroke
            AND NOT has_other_intracranial_hemorrhage
        ) AS subtype_eligible_flag
    FROM labeled
), final_cohort AS (
    SELECT
        'eICU-CRD'::text AS data_source,
        patientunitstayid,
        patienthealthsystemstayid,
        hospitalid,
        hospitaladmitoffset,
        hospitaldischargeoffset,
        hospitaldischargestatus,
        CASE
            WHEN hospitaldischargeoffset >= 1440
                THEN 1 ELSE 0
        END AS landmark_24h_eligible,
        first_qualifying_code_offset,
        first_icd9_stroke_code_offset,
        qualifying_code_in_first_24h,
        primary_icd9_code_in_first_24h,
        audit_phenotype,
        primary_cohort_flag,
        subtype_eligible_flag,
        has_ais,
        has_ich,
        has_sah,
        has_cvst,
        has_unspecified_stroke,
        has_other_intracranial_hemorrhage,
        has_traumatic_intracranial_injury,
        has_ais_icd10_concordant,
        has_ich_icd10_concordant,
        has_sah_icd10_concordant,
        has_unspecified_icd10_concordant,
        has_other_ich_icd10_concordant,
        n_primary_acute_groups,
        CASE
            WHEN lower(coalesce(hospitaldischargestatus, '')) LIKE '%expir%'
                 THEN 1 ELSE 0
        END AS hospital_death_flag,
        CASE
            WHEN lower(coalesce(hospitaldischargestatus, '')) LIKE '%expir%'
                  AND hospitaldischargeoffset >= 0
                  AND hospitaldischargeoffset < 43200
                THEN 1 ELSE 0
        END AS death_30d_inhospital,
        CASE
            WHEN lower(coalesce(hospitaldischargestatus, '')) LIKE '%expir%'
                  AND hospitaldischargeoffset >= 1440
                  AND hospitaldischargeoffset < 43200
                THEN 1 ELSE 0
        END AS death_24h_to_30d_inhospital,
        CASE
            WHEN hospitaldischargestatus IS NULL
                 OR hospitaldischargeoffset IS NULL
                THEN 1 ELSE 0
        END AS outcome_missing_flag
    FROM eligible
)
SELECT *
FROM final_cohort
WHERE primary_cohort_flag
ORDER BY patienthealthsystemstayid, patientunitstayid;
