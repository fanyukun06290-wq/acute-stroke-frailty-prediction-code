, exams AS (
 SELECT s.stay_id,p.physicalexamoffset AS exam_offset,p.physicalexampath AS path
 FROM selected s JOIN eicu_crd.physicalexam p ON p.patientunitstayid=s.stay_id
 JOIN eicu_crd.patient pat ON pat.patientunitstayid=s.stay_id
 WHERE p.physicalexamoffset>=0 AND p.physicalexamoffset<1440
 AND p.physicalexamoffset<pat.unitdischargeoffset
 AND p.physicalexampath LIKE '%/Neurologic/GCS/%'
), parsed AS (
 SELECT stay_id,exam_offset,
 CASE WHEN path ~ '/Eyes Score/[0-9]+$' THEN 'eye'
 WHEN path ~ '/Verbal Score/[0-9]+$' THEN 'verbal'
 WHEN path ~ '/Motor Score/[0-9]+$' THEN 'motor'
 WHEN path ~ '/GCS/[0-9]+$' THEN 'direct_total' ELSE 'quality' END AS component,
 CASE WHEN path ~ '/(Eyes Score|Verbal Score|Motor Score|GCS)/[0-9]+$'
 THEN substring(path FROM '/([0-9]+)$')::float8 ELSE NULL::float8 END AS score,
 (path LIKE '%/Score/estimated due to meds' OR path LIKE '%/Score/unable to score due to meds'
 OR path LIKE '%/Level of Consciousness/sedated') AS explicit_invalid
 FROM exams
)
 SELECT stay_id,CASE WHEN component='quality' THEN 'quality' ELSE 'component' END::text AS kind,
 exam_offset::float8 AS observation_offset,exam_offset::float8 AS available_offset,component,score,
 explicit_invalid,''::text AS airway_type FROM parsed
 UNION ALL
 SELECT s.stay_id,'airway'::text,r.respcarestatusoffset::float8,r.respcarestatusoffset::float8,
 'airway'::text,NULL::float8,FALSE,r.airwaytype
 FROM selected s JOIN eicu_crd.respiratorycare r ON r.patientunitstayid=s.stay_id
 JOIN eicu_crd.patient pat ON pat.patientunitstayid=s.stay_id
 WHERE r.respcarestatusoffset>=0 AND r.respcarestatusoffset<1440
 AND r.respcarestatusoffset<pat.unitdischargeoffset AND nullif(trim(r.airwaytype),'') IS NOT NULL
 ORDER BY stay_id,observation_offset,component
