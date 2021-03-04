context("test-mig_un_fam")

# get data spreadsheet----------------------------------------------------------------
UN_fam <-
  data.frame(
  Age = rep(seq(0,80,5),6),
  Type = c(rep("Family Immigration",17),
           rep("Male Labor Immigration",17),
           rep("Female Labor Immigration",17),
           rep("Family Emigration",17),
           rep("Male Labor Emigration",17),
           rep("Female Labor Emigration",17)),
  Male = c(3082.9,
1339.8,
1027.1,
5237.0,
9426.2,
9059.2,
6861.3,
4758.8,
3195.5,
2120.5,
1401.2,
924.4,
609.6,
401.9,
264.9,
174.6,
115.1,
2136.3,
928.4,
1020.0,
7082.7,
12987.6,
12522.1,
9494.8,
6589.0,
4425.9,
2937.7,
1941.4,
1280.9,
844.7,
556.9,
367.1,
242.0,
159.6,
3945.5,
1714.7,
1029.8,
3513.1,
6102.1,
5827.5,
4403.8,
3050.9,
2047.3,
1358.0,
897.1,
591.7,
390.1,
257.2,
169.5,
111.8,
73.7,
-3082.9,
-1339.8,
-1027.1,
-5237.0,
-9426.2,
-9059.2,
-6861.3,
-4758.8,
-3195.5,
-2120.5,
-1401.2,
-924.4,
-609.6,
-401.9,
-264.9,
-174.6,
-115.1,
-2136.3,
-928.4,
-1020.0,
-7082.7,
-12987.6,
-12522.1,
-9494.8,
-6589.0,
-4425.9,
-2937.7,
-1941.4,
-1280.9,
-844.7,
-556.9,
-367.1,
-242.0,
-159.6,
-3945.5,
-1714.7,
-1029.8,
-3513.1,
-6102.1,
-5827.5,
-4403.8,
-3050.9,
-2047.3,
-1358.0,
-897.1,
-591.7,
-390.1,
-257.2,
-169.5,
-111.8,
-73.7),
Female = c(
  3080.3,
  1338.7,
  1095.0,
  5789.1,
  10073.8,
  9327.4,
  6803.1,
  4543.3,
  2937.5,
  1876.8,
  1194.0,
  758.5,
  481.5,
  305.7,
  194.0,
  123.1,
  78.2,
  2124.3,
  923.2,
  755.2,
  3992.5,
  6947.4,
  6432.7,
  4691.8,
  3133.3,
  2025.8,
  1294.4,
  823.5,
  523.1,
  332.1,
  210.8,
  133.8,
  84.9,
  53.9,
  3974.5,
  1727.3,
  1412.9,
  7469.8,
  12998.4,
  12035.4,
  8778.2,
  5862.3,
  3790.3,
  2421.7,
  1540.7,
  978.7,
  621.3,
  394.4,
  250.3,
  158.9,
  100.9,
  -3080.3,
  -1338.7,
  -1095.0,
  -5789.1,
  -10073.8,
  -9327.4,
  -6803.1,
  -4543.3,
  -2937.5,
  -1876.8,
  -1194.0,
  -758.5,
  -481.5,
  -305.7,
  -194.0,
  -123.1,
  -78.2,
  -2124.3,
  -923.2,
  -755.2,
  -3992.5,
  -6947.4,
  -6432.7,
  -4691.8,
  -3133.3,
  -2025.8,
  -1294.4,
  -823.5,
  -523.1,
  -332.1,
  -210.8,
  -133.8,
  -84.9,
  -53.9,
  -3974.5,
  -1727.3,
  -1412.9,
  -7469.8,
  -12998.4,
  -12035.4,
  -8778.2,
  -5862.3,
  -3790.3,
  -2421.7,
  -1540.7,
  -978.7,
  -621.3,
  -394.4,
  -250.3,
  -158.9,
  -100.9)
) %>%
  rename(Type=2, Age=1) %>%
  as.data.table() %>%
  data.table::melt(id.vars = c("Age","Type"),
                   measure.vars = c("Female","Male"),
                   variable.name = "Sex",
                   value.name = "Prop") %>%
  .[, Prop := Prop / 1e5]


# test --------------------------------------------------------------------

tolerance_admited <- .005
test_that("mig fam works", {
  UN_fam <- UN_fam[Age < 80, ]
  res1 <- mig_un_fam(NM = 1000,  family = "Family", Single = FALSE)$net_migr
  res1 <- res1[res1$age < 80, ]

  expect_equal(
    res1$nm / 1000,
    setDT(UN_fam)[Type == "Family Immigration", .(Prop)][[1]],
    tolerance = tolerance_admited  * 1000
  )

  res2 <- mig_un_fam(NM = -1,  family = "Female Labor", Single = FALSE)$net_migr
  res2 <- res2[res2$age < 80, ]

  expect_equal(
    res2$nm,
    setDT(UN_fam)[Type == "Female Labor Emigration", .(Prop)][[1]],
    tolerance = tolerance_admited
  )

  res3 <- mig_un_fam(NM = -100000,  family = "Male Labor", Single = FALSE)$net_migr
  res3 <- res3[res3$age < 80, ]

  expect_equal(
    res3$nm,
    setDT(UN_fam)[Type == "Male Labor Emigration", .(Prop)][[1]] * 100000,
    tolerance = tolerance_admited * 100000
  )

})

test_that("mig_fam works with OAnew", {

  # Run mig_un_fam with all single ages
  unconstrained <- mig_un_fam(NM = 1000,  family = "Family", OAnew = 120)$net_migr
  # Run with open age group at 100
  constrained <- mig_un_fam(NM = 1000,  family = "Family", OAnew = 100)$net_migr

  # Make sure both calculates sum up to the same when comparing
  # all values above 99
  total_above_100 <- round(sum(unconstrained[age > 99, nm]), 4)
  summed_above_100 <- round(sum(constrained[age == 100, nm]), 4)
  expect_true(total_above_100 == summed_above_100)

})
