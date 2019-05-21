DROP TABLE IF EXISTS overlap;

CREATE TABLE overlap AS
WITH shared_samples AS (
    SELECT DISTINCT sample FROM gdc
    INTERSECT
    SELECT DISTINCT sample FROM disease_wg
    EXCEPT
    -- Exclude samples with duplicated GDC MAFs
    SELECT * FROM (
        VALUES ('C3L-00908'),
               ('C3N-00545'),
               ('C3N-01825')
    )
)
SELECT
    sample,
    chromosome, start_position, end_position, reference_allele, tumor_seq_allele2,
    d.id AS dwg_id, g.id AS gdc_id,
    (CASE WHEN d.id IS NOT NULL AND g.id IS NOT NULL THEN 1 ELSE 0 END) AS shared,
    (CASE WHEN d.id IS NOT NULL AND g.id IS NULL THEN 1 ELSE 0 END) AS only_in_dwg,
    (CASE WHEN d.id IS NULL AND g.id IS NOT NULL THEN 1 ELSE 0 END) AS only_in_gdc
FROM disease_wg d
LEFT JOIN gdc g
    USING (sample, chromosome, start_position, end_position, reference_allele, tumor_seq_allele2)
WHERE sample IN (SELECT sample FROM shared_samples)
UNION ALL
SELECT
    sample,
    chromosome, start_position, end_position, reference_allele, tumor_seq_allele2,
    d.id AS dwg_id, g.id AS gdc_id,
    (CASE WHEN d.id IS NOT NULL AND g.id IS NOT NULL THEN 1 ELSE 0 END) AS shared,
    (CASE WHEN d.id IS NOT NULL AND g.id IS NULL THEN 1 ELSE 0 END) AS only_in_dwg,
    (CASE WHEN d.id IS NULL AND g.id IS NOT NULL THEN 1 ELSE 0 END) AS only_in_gdc
FROM gdc g
LEFT JOIN disease_wg d
    USING (sample, chromosome, start_position, end_position, reference_allele, tumor_seq_allele2)
WHERE d.reference_allele IS NULL
    AND sample IN (SELECT sample FROM shared_samples)
;
