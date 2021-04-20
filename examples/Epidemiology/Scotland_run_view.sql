DROP VIEW IF EXISTS scottish_population_view;
DROP VIEW IF EXISTS pollution_grid_view;

CREATE VIEW scottish_population_view AS
SELECT age_groups
, CAST(substr(age_groups, 4, 5) AS INT) AS age
, (CAST(substr(age_groups, 4, 5) AS INT) / 10) * 10 AS age_aggr
, grid_area
, val
FROM km_age_persons_arr;

CREATE VIEW pollution_grid_view AS
SELECT pollutant
, CAST(substr(grid, 1, instr(grid, '-') - 1) AS INT) AS grid_x
, CAST(substr(grid, instr(grid, '-') + 1) AS INT) AS grid_y
, val
FROM records_pollution_array_arr;
