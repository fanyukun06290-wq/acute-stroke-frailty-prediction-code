/*
ANALYSIS-READY FEATURE DATASET: MIMIC-IV v3.1

Run this file against a PostgreSQL instance containing the MIMIC-IV v3.1 schema.
This is a READ-ONLY query. It creates no database objects and does not
modify any table.

Output: one row per first ICU stay within each hospital admission in the
final primary acute-stroke cohort. The feature policy is fixed as follows:
  - first_6h: first valid value from ICU admission through <6 hours;
  - mean_24h: mean valid value from ICU admission through <24 hours,
              available only for 24-hour landmark-eligible stays;
  - missing_*: explicit missingness indicators for each window.

The 24-hour model must be restricted to landmark_24h_eligible = 1 and use
death_24h_to_30d_inhospital. The 6-hour model may use the full cohort and
death_30d_inhospital.
*/

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
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code IN
                ('43301','43311','43321','43331','43381','43391',
                 '43401','43411','43491'))
            OR (icd_version = 10 AND code LIKE 'I63%'
                AND code NOT LIKE 'I636%')
        ), FALSE) AS has_ais,
        COALESCE(BOOL_OR(
            icd_version = 10 AND code LIKE 'I636%'
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
        ), FALSE) AS has_history_or_sequela,

        /* Chronic-condition proxies used only for frailty-related burden. */
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '250%')
            OR (icd_version = 10 AND substring(code FROM 1 FOR 3)
                IN ('E08','E09','E10','E11','E13'))
        ), FALSE) AS has_diabetes,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code LIKE '428%')
            OR (icd_version = 10 AND code LIKE 'I50%')
        ), FALSE) AS has_heart_failure,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND (code LIKE '491%' OR code LIKE '492%'
                OR code = '496'))
            OR (icd_version = 10 AND code LIKE 'J44%')
        ), FALSE) AS has_copd,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND (code LIKE '290%' OR code LIKE '2941%'
                OR code LIKE '3310%'))
            OR (icd_version = 10 AND (code LIKE 'F01%' OR code LIKE 'F02%'
                OR code LIKE 'F03%' OR code LIKE 'G30%'))
        ), FALSE) AS has_dementia,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND (code LIKE '585%' OR code = '586'))
            OR (icd_version = 10 AND (code LIKE 'N18%' OR code LIKE 'N19%'))
        ), FALSE) AS has_renal_failure,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND code ~ '^(14[0-9]|1[5-9][0-9]|20[0-8])')
            OR (icd_version = 10 AND code LIKE 'C%')
        ), FALSE) AS has_cancer,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND (code LIKE '440%' OR code LIKE '441%'
                OR code LIKE '442%' OR code LIKE '443%'))
            OR (icd_version = 10 AND (code LIKE 'I70%' OR code LIKE 'I71%'
                OR code LIKE 'I72%' OR code LIKE 'I73%' OR code LIKE 'I74%'))
        ), FALSE) AS has_pvd,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND (code LIKE '410%' OR code = '412'))
            OR (icd_version = 10 AND (code LIKE 'I21%' OR code LIKE 'I22%'
                OR code = 'I252'))
        ), FALSE) AS has_prior_mi,
        COALESCE(BOOL_OR(
            (icd_version = 9 AND (code LIKE '438%' OR code = 'V1254'))
            OR (icd_version = 10 AND (code LIKE 'I69%' OR code LIKE 'Z8673%'))
        ), FALSE) AS has_prior_stroke
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
        p.gender,
        (p.anchor_age
            + EXTRACT(YEAR FROM a.admittime)::integer
            - p.anchor_year)::numeric AS age_years,
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
        COALESCE(f.has_tia, FALSE) AS has_tia,
        COALESCE(f.has_history_or_sequela, FALSE)
            AS has_history_or_sequela,
        COALESCE(f.has_diabetes, FALSE) AS has_diabetes,
        COALESCE(f.has_heart_failure, FALSE) AS has_heart_failure,
        COALESCE(f.has_copd, FALSE) AS has_copd,
        COALESCE(f.has_dementia, FALSE) AS has_dementia,
        COALESCE(f.has_renal_failure, FALSE) AS has_renal_failure,
        COALESCE(f.has_cancer, FALSE) AS has_cancer,
        COALESCE(f.has_pvd, FALSE) AS has_pvd,
        COALESCE(f.has_prior_mi, FALSE) AS has_prior_mi,
        COALESCE(f.has_prior_stroke, FALSE) AS has_prior_stroke,
        (f.hadm_id IS NULL)::integer AS diagnosis_record_missing_flag
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
), primary_cohort AS (
    SELECT
        *,
        (n_primary_acute_groups = 1
         AND NOT has_traumatic_intracranial_injury
         AND NOT has_cvst
         AND NOT has_unspecified_stroke
         AND NOT has_other_intracranial_hemorrhage) AS subtype_eligible_flag,
        CASE WHEN dischtime >= icu_intime + INTERVAL '24 hours'
            THEN 1 ELSE 0 END AS landmark_24h_eligible,
        CASE WHEN hospital_expire_flag = 1
                  AND deathtime IS NOT NULL
                  AND deathtime >= icu_intime
                  AND deathtime < icu_intime + INTERVAL '30 days'
            THEN 1 ELSE 0 END AS death_30d_inhospital,
        CASE WHEN hospital_expire_flag = 1
                  AND deathtime IS NOT NULL
                  AND deathtime >= icu_intime + INTERVAL '24 hours'
                  AND deathtime < icu_intime + INTERVAL '30 days'
            THEN 1 ELSE 0 END AS death_24h_to_30d_inhospital,
        CASE WHEN hospital_expire_flag IS NULL
                  OR (hospital_expire_flag = 1 AND deathtime IS NULL)
            THEN 1 ELSE 0 END AS outcome_missing_flag
    FROM labeled
    WHERE n_primary_acute_groups > 0
      AND NOT has_traumatic_intracranial_injury
), lab_map(variable, itemid) AS (
    VALUES
        ('glucose', 50809), ('glucose', 50931), ('glucose', 52569),
        ('sodium', 50983), ('sodium', 52623),
        ('potassium', 50971), ('potassium', 52610),
        ('chloride', 50902), ('chloride', 52535),
        ('bicarbonate', 50882),
        ('bun', 51006), ('bun', 52647),
        ('creatinine', 50912), ('creatinine', 52546),
        ('calcium_total', 50893),
        ('wbc', 51300), ('wbc', 51301),
        ('platelets', 51265),
        ('hemoglobin', 51222),
        ('hematocrit', 51221)
), lab_bounds(variable, lower_bound, upper_bound) AS (
    -- Outcome-blind broad screening bounds; verify units/item meanings before final use.
    VALUES
        ('glucose', 20, 1000), ('sodium', 100, 200),
        ('potassium', 1.5, 10), ('chloride', 50, 150),
        ('bicarbonate', 5, 60), ('bun', 1, 300),
        ('creatinine', 0.1, 30), ('calcium_total', 4, 15),
        ('wbc', 0.1, 300), ('platelets', 10, 2000),
        ('hemoglobin', 2, 25), ('hematocrit', 5, 75)
), lab_long AS (
    SELECT
        c.stay_id,
        lm.variable,
        le.charttime AS measured_at,
        EXTRACT(EPOCH FROM (le.charttime - c.icu_intime)) / 60.0
            AS offset_min,
        le.valuenum::numeric AS value
    FROM primary_cohort AS c
    INNER JOIN mimiciv_hosp.labevents AS le
        ON le.hadm_id = c.hadm_id
       AND le.charttime >= c.icu_intime
       AND le.charttime < c.icu_intime + INTERVAL '24 hours'
       AND le.valuenum IS NOT NULL
    INNER JOIN lab_map AS lm
        ON lm.itemid = le.itemid
    INNER JOIN lab_bounds AS lb
        ON lb.variable = lm.variable
       AND le.valuenum::numeric BETWEEN lb.lower_bound AND lb.upper_bound
), lab_summary AS (
    SELECT
        c.stay_id,
        lm.variable,
        (array_agg(l.value ORDER BY l.measured_at)
            FILTER (WHERE l.offset_min < 360))[1] AS first_6h,
        CASE WHEN c.landmark_24h_eligible = 1
            THEN AVG(l.value) FILTER (WHERE l.offset_min < 1440)
        END AS mean_24h,
        CASE WHEN COUNT(l.value) FILTER (WHERE l.offset_min < 360) > 0
            THEN 0 ELSE 1 END AS missing_6h,
        CASE WHEN c.landmark_24h_eligible = 1 THEN
            CASE WHEN COUNT(l.value) FILTER (WHERE l.offset_min < 1440) > 0
                THEN 0 ELSE 1 END
        END AS missing_24h
    FROM primary_cohort AS c
    CROSS JOIN (SELECT DISTINCT variable FROM lab_map) AS lm
    LEFT JOIN lab_long AS l
        ON l.stay_id = c.stay_id
       AND l.variable = lm.variable
    GROUP BY c.stay_id, c.landmark_24h_eligible, lm.variable
), lab_wide AS (
    SELECT
        stay_id AS lab_stay_id,
        MAX(first_6h) FILTER (WHERE variable = 'glucose')
            AS glucose_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'glucose')
            AS glucose_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'glucose')
            AS glucose_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'glucose')
            AS glucose_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'sodium') AS sodium_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'sodium') AS sodium_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'sodium') AS sodium_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'sodium') AS sodium_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'potassium') AS potassium_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'potassium') AS potassium_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'potassium') AS potassium_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'potassium') AS potassium_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'chloride') AS chloride_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'chloride') AS chloride_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'chloride') AS chloride_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'chloride') AS chloride_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'bicarbonate') AS bicarbonate_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'bicarbonate') AS bicarbonate_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'bicarbonate') AS bicarbonate_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'bicarbonate') AS bicarbonate_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'bun') AS bun_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'bun') AS bun_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'bun') AS bun_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'bun') AS bun_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'creatinine') AS creatinine_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'creatinine') AS creatinine_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'creatinine') AS creatinine_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'creatinine') AS creatinine_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'calcium_total') AS calcium_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'calcium_total') AS calcium_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'calcium_total') AS calcium_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'calcium_total') AS calcium_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'wbc') AS wbc_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'wbc') AS wbc_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'wbc') AS wbc_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'wbc') AS wbc_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'platelets') AS platelets_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'platelets') AS platelets_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'platelets') AS platelets_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'platelets') AS platelets_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'hemoglobin') AS hemoglobin_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'hemoglobin') AS hemoglobin_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'hemoglobin') AS hemoglobin_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'hemoglobin') AS hemoglobin_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'hematocrit') AS hematocrit_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'hematocrit') AS hematocrit_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'hematocrit') AS hematocrit_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'hematocrit') AS hematocrit_missing_24h
    FROM lab_summary
    GROUP BY stay_id
), vital_map(variable, itemid) AS (
    VALUES
        ('heart_rate', 220045),
        ('respiratory_rate', 220210),
        ('spo2', 220277),
        ('sbp', 220050), ('sbp', 220179),
        ('dbp', 220051), ('dbp', 220180),
        ('map', 220052), ('map', 220181)
), vital_bounds(variable, lower_bound, upper_bound) AS (
    -- Outcome-blind broad screening bounds; verify item meanings before final use.
    VALUES
        ('heart_rate', 20, 250), ('respiratory_rate', 4, 100),
        ('spo2', 50, 100), ('sbp', 50, 300),
        ('dbp', 20, 200), ('map', 30, 250)
), vital_long AS (
    SELECT
        c.stay_id,
        vm.variable,
        ce.charttime AS measured_at,
        EXTRACT(EPOCH FROM (ce.charttime - c.icu_intime)) / 60.0
            AS offset_min,
        ce.valuenum::numeric AS value
    FROM primary_cohort AS c
    INNER JOIN mimiciv_icu.chartevents AS ce
        ON ce.stay_id = c.stay_id
       AND ce.charttime >= c.icu_intime
       AND ce.charttime < c.icu_intime + INTERVAL '24 hours'
       AND ce.valuenum IS NOT NULL
    INNER JOIN vital_map AS vm
        ON vm.itemid = ce.itemid
    INNER JOIN vital_bounds AS vb
        ON vb.variable = vm.variable
       AND ce.valuenum::numeric BETWEEN vb.lower_bound AND vb.upper_bound
), vital_summary AS (
    SELECT
        c.stay_id,
        vm.variable,
        (array_agg(v.value ORDER BY v.measured_at)
            FILTER (WHERE v.offset_min < 360))[1] AS first_6h,
        CASE WHEN c.landmark_24h_eligible = 1
            THEN AVG(v.value) FILTER (WHERE v.offset_min < 1440)
        END AS mean_24h,
        CASE WHEN COUNT(v.value) FILTER (WHERE v.offset_min < 360) > 0
            THEN 0 ELSE 1 END AS missing_6h,
        CASE WHEN c.landmark_24h_eligible = 1 THEN
            CASE WHEN COUNT(v.value) FILTER (WHERE v.offset_min < 1440) > 0
                THEN 0 ELSE 1 END
        END AS missing_24h
    FROM primary_cohort AS c
    CROSS JOIN (SELECT DISTINCT variable FROM vital_map) AS vm
    LEFT JOIN vital_long AS v
        ON v.stay_id = c.stay_id
       AND v.variable = vm.variable
    GROUP BY c.stay_id, c.landmark_24h_eligible, vm.variable
), vital_wide AS (
    SELECT
        stay_id AS vital_stay_id,
        MAX(first_6h) FILTER (WHERE variable = 'heart_rate') AS heart_rate_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'heart_rate') AS heart_rate_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'heart_rate') AS heart_rate_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'heart_rate') AS heart_rate_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'respiratory_rate') AS respiratory_rate_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'respiratory_rate') AS respiratory_rate_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'respiratory_rate') AS respiratory_rate_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'respiratory_rate') AS respiratory_rate_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'spo2') AS spo2_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'spo2') AS spo2_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'spo2') AS spo2_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'spo2') AS spo2_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'sbp') AS sbp_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'sbp') AS sbp_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'sbp') AS sbp_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'sbp') AS sbp_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'dbp') AS dbp_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'dbp') AS dbp_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'dbp') AS dbp_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'dbp') AS dbp_missing_24h,
        MAX(first_6h) FILTER (WHERE variable = 'map') AS map_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'map') AS map_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'map') AS map_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'map') AS map_missing_24h
    FROM vital_summary
    GROUP BY stay_id
)
SELECT
    'MIMIC-IV'::text AS data_source,
    c.subject_id,
    c.hadm_id,
    c.stay_id,
    c.age_years,
    upper(trim(c.gender)) AS sex,
    c.icu_intime,
    c.landmark_24h_eligible,
    c.audit_phenotype,
    c.subtype_eligible_flag,
    c.has_ais,
    c.has_ich,
    c.has_sah,
    c.has_cvst,
    c.has_unspecified_stroke,
    c.has_other_intracranial_hemorrhage,
    c.n_primary_acute_groups,
    c.has_diabetes,
    c.has_heart_failure,
    c.has_copd,
    c.has_dementia,
    c.has_renal_failure,
    c.has_cancer,
    c.has_pvd,
    c.has_prior_mi,
    c.has_prior_stroke,
    (c.has_diabetes::integer + c.has_heart_failure::integer
     + c.has_copd::integer + c.has_dementia::integer
     + c.has_renal_failure::integer + c.has_cancer::integer
     + c.has_pvd::integer + c.has_prior_mi::integer
     + c.has_prior_stroke::integer) AS frailty_comorbidity_count,
    c.diagnosis_record_missing_flag,
    c.death_30d_inhospital,
    c.death_24h_to_30d_inhospital,
    c.outcome_missing_flag,
    lw.*,
    vw.*
FROM primary_cohort AS c
LEFT JOIN lab_wide AS lw ON lw.lab_stay_id = c.stay_id
LEFT JOIN vital_wide AS vw ON vw.vital_stay_id = c.stay_id
ORDER BY c.hadm_id, c.stay_id;
