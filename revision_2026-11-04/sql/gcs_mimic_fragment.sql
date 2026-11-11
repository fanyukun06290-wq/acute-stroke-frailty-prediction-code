
 SELECT c.stay_id,'component'::text AS kind,
 extract(epoch FROM(c.charttime-i.intime))/60.0 AS observation_offset,
 extract(epoch FROM(c.storetime-i.intime))/60.0 AS available_offset,
 CASE c.itemid WHEN 220739 THEN 'eye' WHEN 223900 THEN 'verbal' ELSE 'motor' END AS component,
 c.valuenum::float8 AS score,
 coalesce(c.value,'') ~* '(ETT|ET/|Trach|unable|unassess|not.?test|non.?test|sedat)' AS explicit_invalid,
 ''::text AS airway_type
 FROM selected s JOIN mimiciv_icu.icustays i USING(stay_id)
 JOIN mimiciv_icu.chartevents c ON c.stay_id=s.stay_id
 WHERE c.itemid IN(220739,223900,223901) AND c.charttime>=i.intime
 AND c.charttime<i.intime+interval '24 hours' AND c.charttime<i.outtime
 ORDER BY c.stay_id,c.charttime,c.itemid,c.storetime
