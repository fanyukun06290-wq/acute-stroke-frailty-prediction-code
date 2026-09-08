/*
ANALYSIS-READY FEATURE DATASET: eICU-CRD v2.0

Run this file against a PostgreSQL instance containing the eICU-CRD v2.0 schema.
This is a READ-ONLY query. It creates no database objects and does not
modify any table.

Output: one row per first ICU unit within each hospital stay in the final
primary acute-stroke cohort. The feature policy is fixed as follows:
  - first_6h: first valid value from ICU admission through <6 hours;
  - mean_24h: mean valid value from ICU admission through <24 hours,
              available only for 24-hour landmark-eligible stays;
  - missing_*: explicit missingness indicators for each window.

The 24-hour model must be restricted to landmark_24h_eligible = 1 and use
death_24h_to_30d_inhospital. The 6-hour model may use the full cohort and
death_30d_inhospital.
*/

WITH first_unit AS (
    SELECT
        p.patientunitstayid,
        p.patienthealthsystemstayid,
        p.hospitalid,
        p.age,
        p.gender,
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
), past_history_flags AS (
    SELECT
        ph.patientunitstayid,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%diabetes%'
            AND lower(ph.pasthistorypath) NOT LIKE '%none%'
        ), FALSE) AS has_diabetes,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%congestive heart failure%'
            AND lower(ph.pasthistorypath) NOT LIKE '%none%'
        ), FALSE) AS has_heart_failure,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%/copd/%'
            AND lower(ph.pasthistorypath) NOT LIKE '%none%'
        ), FALSE) AS has_copd,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%/dementia/%'
            AND lower(ph.pasthistorypath) NOT LIKE '%none%'
        ), FALSE) AS has_dementia,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%/renal failure/%'
            OR lower(ph.pasthistorypath) LIKE '%/renal insufficiency/%'
        ), FALSE) AS has_renal_failure,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%/cancer/%'
            AND lower(ph.pasthistorypath) NOT LIKE '%cancer therapy%'
            AND lower(ph.pasthistorypath) NOT LIKE '%/none%'
        ), FALSE) AS has_cancer,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%peripheral vascular disease%'
            AND lower(ph.pasthistorypath) NOT LIKE '%none%'
        ), FALSE) AS has_pvd,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%myocardial infarction%'
            AND lower(ph.pasthistorypath) NOT LIKE '%none%'
        ), FALSE) AS has_prior_mi,
        COALESCE(BOOL_OR(
            lower(ph.pasthistorypath) LIKE '%/strokes/%'
            AND lower(ph.pasthistorypath) NOT LIKE '%none%'
        ), FALSE) AS has_prior_stroke
    FROM eicu_crd.pasthistory AS ph
    WHERE ph.pasthistoryoffset >= 0
      AND ph.pasthistoryoffset < 1440
    GROUP BY ph.patientunitstayid
), cohort_base AS (
    SELECT
        p.patientunitstayid,
        p.patienthealthsystemstayid,
        p.hospitalid,
        p.hospitaladmitoffset,
        p.hospitaldischargeoffset,
        p.hospitaldischargestatus,
        CASE
            WHEN trim(p.age) ~ '^[0-9]+$' THEN trim(p.age)::numeric
            WHEN trim(p.age) LIKE '>%' THEN 90::numeric
            ELSE NULL::numeric
        END AS age_years,
        CASE
            WHEN upper(trim(p.gender)) IN ('MALE','M') THEN 'M'
            WHEN upper(trim(p.gender)) IN ('FEMALE','F') THEN 'F'
            ELSE NULL
        END AS sex,
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
        f.first_qualifying_code_offset,
        COALESCE(ph.has_diabetes, FALSE) AS has_diabetes,
        COALESCE(ph.has_heart_failure, FALSE) AS has_heart_failure,
        COALESCE(ph.has_copd, FALSE) AS has_copd,
        COALESCE(ph.has_dementia, FALSE) AS has_dementia,
        COALESCE(ph.has_renal_failure, FALSE) AS has_renal_failure,
        COALESCE(ph.has_cancer, FALSE) AS has_cancer,
        COALESCE(ph.has_pvd, FALSE) AS has_pvd,
        COALESCE(ph.has_prior_mi, FALSE) AS has_prior_mi,
        COALESCE(ph.has_prior_stroke, FALSE) AS has_prior_stroke,
        (ph.patientunitstayid IS NULL)::integer
            AS past_history_record_missing_flag
    FROM first_unit AS p
    LEFT JOIN flags AS f
        ON f.patientunitstayid = p.patientunitstayid
    LEFT JOIN past_history_flags AS ph
        ON ph.patientunitstayid = p.patientunitstayid
    WHERE p.unit_order = 1
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
        CASE WHEN hospitaldischargeoffset >= 1440
            THEN 1 ELSE 0 END AS landmark_24h_eligible,
        CASE WHEN lower(coalesce(hospitaldischargestatus, '')) LIKE '%expir%'
                  AND hospitaldischargeoffset >= 0
                  AND hospitaldischargeoffset < 43200
            THEN 1 ELSE 0 END AS death_30d_inhospital,
        CASE WHEN lower(coalesce(hospitaldischargestatus, '')) LIKE '%expir%'
                  AND hospitaldischargeoffset >= 1440
                  AND hospitaldischargeoffset < 43200
            THEN 1 ELSE 0 END AS death_24h_to_30d_inhospital,
        CASE WHEN hospitaldischargestatus IS NULL
                  OR hospitaldischargeoffset IS NULL
            THEN 1 ELSE 0 END AS outcome_missing_flag
    FROM labeled
    WHERE n_primary_acute_groups > 0
      AND NOT has_traumatic_intracranial_injury
), lab_map(variable, labname) AS (
    VALUES
        ('glucose', 'glucose'), ('glucose', 'bedside glucose'),
        ('sodium', 'sodium'),
        ('potassium', 'potassium'),
        ('chloride', 'chloride'),
        ('bicarbonate', 'bicarbonate'),
        ('bun', 'bun'),
        ('creatinine', 'creatinine'),
        ('calcium_total', 'calcium'),
        ('wbc', 'wbc x 1000'),
        ('platelets', 'platelets x 1000'),
        ('hemoglobin', 'hgb'),
        ('hematocrit', 'hct')
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
        c.patientunitstayid,
        lm.variable,
        l.labresultoffset::numeric AS measured_at,
        l.labresult::numeric AS value
    FROM primary_cohort AS c
    INNER JOIN eicu_crd.lab AS l
        ON l.patientunitstayid = c.patientunitstayid
       AND l.labresultoffset >= 0
       AND l.labresultoffset < 1440
       AND l.labresult IS NOT NULL
    INNER JOIN lab_map AS lm
        ON trim(lower(l.labname)) = lm.labname
    INNER JOIN lab_bounds AS lb
        ON lb.variable = lm.variable
       AND l.labresult::numeric BETWEEN lb.lower_bound AND lb.upper_bound
), lab_summary AS (
    SELECT
        c.patientunitstayid,
        lm.variable,
        (array_agg(l.value ORDER BY l.measured_at)
            FILTER (WHERE l.measured_at < 360))[1] AS first_6h,
        CASE WHEN c.landmark_24h_eligible = 1
            THEN AVG(l.value) FILTER (WHERE l.measured_at < 1440)
        END AS mean_24h,
        CASE WHEN COUNT(l.value) FILTER (WHERE l.measured_at < 360) > 0
            THEN 0 ELSE 1 END AS missing_6h,
        CASE WHEN c.landmark_24h_eligible = 1 THEN
            CASE WHEN COUNT(l.value) FILTER (WHERE l.measured_at < 1440) > 0
                THEN 0 ELSE 1 END
        END AS missing_24h
    FROM primary_cohort AS c
    CROSS JOIN (SELECT DISTINCT variable FROM lab_map) AS lm
    LEFT JOIN lab_long AS l
        ON l.patientunitstayid = c.patientunitstayid
       AND l.variable = lm.variable
    GROUP BY c.patientunitstayid, c.landmark_24h_eligible, lm.variable
), lab_wide AS (
    SELECT
        patientunitstayid AS lab_patientunitstayid,
        MAX(first_6h) FILTER (WHERE variable = 'glucose') AS glucose_first_6h,
        MAX(mean_24h) FILTER (WHERE variable = 'glucose') AS glucose_mean_24h,
        MAX(missing_6h) FILTER (WHERE variable = 'glucose') AS glucose_missing_6h,
        MAX(missing_24h) FILTER (WHERE variable = 'glucose') AS glucose_missing_24h,
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
    GROUP BY patientunitstayid
), vital_bounds(variable, lower_bound, upper_bound) AS (
    -- Outcome-blind broad screening bounds; verify units/item meanings before final use.
    VALUES
        ('heart_rate', 20, 250), ('respiratory_rate', 4, 100),
        ('spo2', 50, 100), ('sbp', 50, 300),
        ('dbp', 20, 200), ('map', 30, 250)
), vital_long AS (
    SELECT
        c.patientunitstayid,
        v.variable,
        vp.observationoffset::numeric AS measured_at,
        v.value::numeric AS value
    FROM primary_cohort AS c
    INNER JOIN eicu_crd.vitalperiodic AS vp
        ON vp.patientunitstayid = c.patientunitstayid
       AND vp.observationoffset >= 0
       AND vp.observationoffset < 1440
    CROSS JOIN LATERAL (
        VALUES
            ('heart_rate', vp.heartrate::numeric),
            ('respiratory_rate', vp.respiration::numeric),
            ('spo2', vp.sao2::numeric),
            ('sbp', vp.systemicsystolic::numeric),
            ('dbp', vp.systemicdiastolic::numeric),
            ('map', vp.systemicmean::numeric)
    ) AS v(variable, value)
    INNER JOIN vital_bounds AS vb ON vb.variable = v.variable
    WHERE v.value IS NOT NULL
      AND v.value BETWEEN vb.lower_bound AND vb.upper_bound
    UNION ALL
    SELECT
        c.patientunitstayid,
        v.variable,
        va.observationoffset::numeric AS measured_at,
        v.value::numeric AS value
    FROM primary_cohort AS c
    INNER JOIN eicu_crd.vitalaperiodic AS va
        ON va.patientunitstayid = c.patientunitstayid
       AND va.observationoffset >= 0
       AND va.observationoffset < 1440
    CROSS JOIN LATERAL (
        VALUES
            ('sbp', va.noninvasivesystolic::numeric),
            ('dbp', va.noninvasivediastolic::numeric),
            ('map', va.noninvasivemean::numeric)
    ) AS v(variable, value)
    INNER JOIN vital_bounds AS vb ON vb.variable = v.variable
    WHERE v.value IS NOT NULL
      AND v.value BETWEEN vb.lower_bound AND vb.upper_bound
), vital_summary AS (
    SELECT
        c.patientunitstayid,
        vm.variable,
        (array_agg(v.value ORDER BY v.measured_at)
            FILTER (WHERE v.measured_at < 360))[1] AS first_6h,
        CASE WHEN c.landmark_24h_eligible = 1
            THEN AVG(v.value) FILTER (WHERE v.measured_at < 1440)
        END AS mean_24h,
        CASE WHEN COUNT(v.value) FILTER (WHERE v.measured_at < 360) > 0
            THEN 0 ELSE 1 END AS missing_6h,
        CASE WHEN c.landmark_24h_eligible = 1 THEN
            CASE WHEN COUNT(v.value) FILTER (WHERE v.measured_at < 1440) > 0
                THEN 0 ELSE 1 END
        END AS missing_24h
    FROM primary_cohort AS c
    CROSS JOIN (VALUES
        ('heart_rate'), ('respiratory_rate'), ('spo2'),
        ('sbp'), ('dbp'), ('map')
    ) AS vm(variable)
    LEFT JOIN vital_long AS v
        ON v.patientunitstayid = c.patientunitstayid
       AND v.variable = vm.variable
    GROUP BY c.patientunitstayid, c.landmark_24h_eligible, vm.variable
), vital_wide AS (
    SELECT
        patientunitstayid AS vital_patientunitstayid,
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
    GROUP BY patientunitstayid
)
SELECT
    'eICU-CRD'::text AS data_source,
    c.patientunitstayid,
    c.patienthealthsystemstayid,
    c.hospitalid,
    c.age_years,
    c.sex,
    c.hospitaladmitoffset,
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
    c.past_history_record_missing_flag,
    c.death_30d_inhospital,
    c.death_24h_to_30d_inhospital,
    c.outcome_missing_flag,
    lw.*,
    vw.*
FROM primary_cohort AS c
LEFT JOIN lab_wide AS lw
    ON lw.lab_patientunitstayid = c.patientunitstayid
LEFT JOIN vital_wide AS vw
    ON vw.vital_patientunitstayid = c.patientunitstayid
ORDER BY c.patienthealthsystemstayid, c.patientunitstayid;
