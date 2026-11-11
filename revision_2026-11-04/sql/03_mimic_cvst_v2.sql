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
), feature_cohort AS (
 SELECT * FROM primary_cohort WHERE landmark_6h_eligible=1 OR landmark_24h_eligible=1
), bounds(variable,lower_bound,upper_bound) AS (VALUES ('glucose',20,1000),
('sodium',100,200),
('potassium',1.5,10),
('chloride',50,150),
('bicarbonate',5,60),
('bun',1,300),
('creatinine',0.1,30),
('calcium',4,15),
('wbc',0.1,300),
('platelets',10,2000),
('hemoglobin',2,25),
('hematocrit',5,75),
('heart_rate',20,250),
('respiratory_rate',4,100),
('spo2',50,100),
('sbp',50,300),
('dbp',20,200),
('map',30,250)), lab_map(variable,itemid) AS (
    VALUES
        ('glucose', 50809), ('glucose', 50931), ('glucose', 52569),
        ('sodium', 50983), ('sodium', 52623),
        ('potassium', 50971), ('potassium', 52610),
        ('chloride', 50902), ('chloride', 52535),
        ('bicarbonate', 50882),
        ('bun', 51006), ('bun', 52647),
        ('creatinine', 50912), ('creatinine', 52546),
        ('calcium', 50893),
        ('wbc', 51300), ('wbc', 51301),
        ('platelets', 51265),
        ('hemoglobin', 51222),
        ('hematocrit', 51221)
), vital_map(variable,itemid) AS (
    VALUES
        ('heart_rate', 220045),
        ('respiratory_rate', 220210),
        ('spo2', 220277),
        ('sbp', 220050), ('sbp', 220179),
        ('dbp', 220051), ('dbp', 220180),
        ('map', 220052), ('map', 220181)
), raw_long AS (
 SELECT c.stay_id,lm.variable,
 extract(epoch FROM (le.charttime-c.icu_intime))/60.0 AS measured_min,
 CASE WHEN le.storetime IS NOT NULL THEN
 extract(epoch FROM (GREATEST(le.charttime,le.storetime)-c.icu_intime))/60.0 END AS available_min,
 le.valuenum::numeric AS value,le.labevent_id::text AS event_key,
 le.valueuom::text AS unit,'lab'::text AS event_source
 FROM feature_cohort c JOIN mimiciv_hosp.labevents le ON le.hadm_id=c.hadm_id
 AND le.charttime>=c.icu_intime AND le.charttime<c.icu_intime+interval '24 hours'
 AND le.charttime<c.effective_icu_end_at
 JOIN lab_map lm ON lm.itemid=le.itemid WHERE le.valuenum IS NOT NULL
 UNION ALL
 SELECT c.stay_id,vm.variable,
 extract(epoch FROM (ce.charttime-c.icu_intime))/60.0,
 CASE WHEN ce.storetime IS NOT NULL THEN
 extract(epoch FROM (GREATEST(ce.charttime,ce.storetime)-c.icu_intime))/60.0 END,
 ce.valuenum::numeric,ce.itemid::text,ce.valueuom::text,'chart_vital'::text
 FROM feature_cohort c JOIN mimiciv_icu.chartevents ce ON ce.stay_id=c.stay_id
 AND ce.charttime>=c.icu_intime AND ce.charttime<c.icu_intime+interval '24 hours'
 AND ce.charttime<c.effective_icu_end_at
 JOIN vital_map vm ON vm.itemid=ce.itemid WHERE ce.valuenum IS NOT NULL
), valid_long AS (
 SELECT r.* FROM raw_long r JOIN bounds b USING(variable)
 WHERE r.value BETWEEN b.lower_bound AND b.upper_bound
), features AS (
 SELECT c.stay_id AS feature_stay_id,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='glucose' AND l.available_min<360 AND l.measured_min<360))[1] END AS glucose_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='glucose' AND l.available_min<1440 AND l.measured_min<1440) END AS glucose_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='glucose' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS glucose_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='glucose' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS glucose_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='sodium' AND l.available_min<360 AND l.measured_min<360))[1] END AS sodium_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='sodium' AND l.available_min<1440 AND l.measured_min<1440) END AS sodium_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='sodium' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS sodium_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='sodium' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS sodium_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='potassium' AND l.available_min<360 AND l.measured_min<360))[1] END AS potassium_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='potassium' AND l.available_min<1440 AND l.measured_min<1440) END AS potassium_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='potassium' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS potassium_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='potassium' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS potassium_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='chloride' AND l.available_min<360 AND l.measured_min<360))[1] END AS chloride_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='chloride' AND l.available_min<1440 AND l.measured_min<1440) END AS chloride_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='chloride' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS chloride_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='chloride' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS chloride_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='bicarbonate' AND l.available_min<360 AND l.measured_min<360))[1] END AS bicarbonate_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='bicarbonate' AND l.available_min<1440 AND l.measured_min<1440) END AS bicarbonate_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='bicarbonate' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS bicarbonate_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='bicarbonate' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS bicarbonate_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='bun' AND l.available_min<360 AND l.measured_min<360))[1] END AS bun_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='bun' AND l.available_min<1440 AND l.measured_min<1440) END AS bun_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='bun' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS bun_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='bun' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS bun_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='creatinine' AND l.available_min<360 AND l.measured_min<360))[1] END AS creatinine_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='creatinine' AND l.available_min<1440 AND l.measured_min<1440) END AS creatinine_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='creatinine' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS creatinine_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='creatinine' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS creatinine_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='calcium' AND l.available_min<360 AND l.measured_min<360))[1] END AS calcium_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='calcium' AND l.available_min<1440 AND l.measured_min<1440) END AS calcium_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='calcium' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS calcium_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='calcium' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS calcium_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='wbc' AND l.available_min<360 AND l.measured_min<360))[1] END AS wbc_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='wbc' AND l.available_min<1440 AND l.measured_min<1440) END AS wbc_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='wbc' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS wbc_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='wbc' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS wbc_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='platelets' AND l.available_min<360 AND l.measured_min<360))[1] END AS platelets_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='platelets' AND l.available_min<1440 AND l.measured_min<1440) END AS platelets_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='platelets' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS platelets_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='platelets' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS platelets_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='hemoglobin' AND l.available_min<360 AND l.measured_min<360))[1] END AS hemoglobin_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='hemoglobin' AND l.available_min<1440 AND l.measured_min<1440) END AS hemoglobin_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='hemoglobin' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS hemoglobin_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='hemoglobin' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS hemoglobin_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='hematocrit' AND l.available_min<360 AND l.measured_min<360))[1] END AS hematocrit_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='hematocrit' AND l.available_min<1440 AND l.measured_min<1440) END AS hematocrit_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='hematocrit' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS hematocrit_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='hematocrit' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS hematocrit_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='heart_rate' AND l.available_min<360 AND l.measured_min<360))[1] END AS heart_rate_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='heart_rate' AND l.available_min<1440 AND l.measured_min<1440) END AS heart_rate_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='heart_rate' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS heart_rate_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='heart_rate' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS heart_rate_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='respiratory_rate' AND l.available_min<360 AND l.measured_min<360))[1] END AS respiratory_rate_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='respiratory_rate' AND l.available_min<1440 AND l.measured_min<1440) END AS respiratory_rate_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='respiratory_rate' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS respiratory_rate_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='respiratory_rate' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS respiratory_rate_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='spo2' AND l.available_min<360 AND l.measured_min<360))[1] END AS spo2_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='spo2' AND l.available_min<1440 AND l.measured_min<1440) END AS spo2_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='spo2' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS spo2_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='spo2' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS spo2_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='sbp' AND l.available_min<360 AND l.measured_min<360))[1] END AS sbp_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='sbp' AND l.available_min<1440 AND l.measured_min<1440) END AS sbp_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='sbp' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS sbp_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='sbp' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS sbp_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='dbp' AND l.available_min<360 AND l.measured_min<360))[1] END AS dbp_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='dbp' AND l.available_min<1440 AND l.measured_min<1440) END AS dbp_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='dbp' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS dbp_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='dbp' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS dbp_missing_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (array_agg(l.value ORDER BY l.measured_min,l.available_min,l.event_key,l.value) FILTER(WHERE l.variable='map' AND l.available_min<360 AND l.measured_min<360))[1] END AS map_first_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN avg(l.value) FILTER(WHERE l.variable='map' AND l.available_min<1440 AND l.measured_min<1440) END AS map_mean_24h,
CASE WHEN c.landmark_6h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='map' AND l.available_min<360 AND l.measured_min<360)=0)::integer END AS map_missing_6h,
CASE WHEN c.landmark_24h_eligible=1 THEN (count(l.value) FILTER(WHERE l.variable='map' AND l.available_min<1440 AND l.measured_min<1440)=0)::integer END AS map_missing_24h,
count(l.value) FILTER(WHERE l.available_min IS NULL) AS audit_rows_unknown_recording_time,
count(l.value) FILTER(WHERE l.event_source='lab' AND l.measured_min<360 AND l.available_min>=360) AS audit_lab_rows_recorded_after6,
count(l.value) FILTER(WHERE l.event_source='lab' AND l.measured_min<1440 AND l.available_min>=1440) AS audit_lab_rows_recorded_after24
 FROM feature_cohort c LEFT JOIN valid_long l USING(stay_id)
 GROUP BY c.stay_id,c.landmark_6h_eligible,c.landmark_24h_eligible
)
SELECT c.*,f.* FROM primary_cohort c LEFT JOIN features f ON f.feature_stay_id=c.stay_id
ORDER BY c.stay_id;
