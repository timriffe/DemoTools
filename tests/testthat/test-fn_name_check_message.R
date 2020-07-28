# context("warnings based on name of called function")
# 
# test_that("calling 'Whipple' produces a warning", {
# 
#     expect_warning(Whipple(Value = 25:65, Age = 25:65), "please use")
#     expect_warning(Whipple(Value = 25:65, Age = 25:65), "please use")
#     expect_warning(do.call("Whipple",
#                            args = list(Value = 25:65, Age = 25:65)), "please use")
# 
#     ## Invoking 'do.call' in these forms does not let 'match.call' see
#     ## the function name so any use of the old function name will not
#     ## be caught.
#     expect_warning(do.call(match.fun("Whipple"),
#                            args = list(Value = 25:65, Age = 25:65)), NA)
#     expect_warning(do.call(Whipple,
#                            args = list(Value = 25:65, Age = 25:65)), NA)
#     expect_warning(do.call(match.fun(Whipple),
#                            args = list(Value = 25:65, Age = 25:65)), NA)
# 
# })
# 
# test_that("calling 'DemoTools::Whipple' proudces a warning", {
# 
#     expect_warning(DemoTools::Whipple(Value = 25:65, Age = 25:65), "please use")
#     expect_warning(DemoTools::Whipple(Value = 25:65, Age = 25:65), "please use")
# 
#     ## Invoking 'do.call' in these forms does not let 'match.call' see
#     ## the function name so any use of the old function name will not
#     ## be caught.
#     expect_warning(do.call(DemoTools::Whipple,
#                            args = list(Value = 25:65, Age = 25:65)), NA)
#     expect_warning(do.call(match.fun(DemoTools::Whipple),
#                            args = list(Value = 25:65, Age = 25:65)), NA)
# })
# 
# 
# test_that("calling 'check_heaping_whipple' does not produce a warning", {
# 
#     expect_warning(check_heaping_whipple(Value = 25:65, Age = 25:65), NA)
#     expect_warning(check_heaping_whipple(Value = 25:65, Age = 25:65), NA)
#     expect_warning(do.call("check_heaping_whipple",
#                            args = list(Value = 25:65, Age = 25:65)), NA)
# 
#     ## Invoking 'do.call' in these forms does not let 'match.call' see
#     ## the function name so any use of the old function name will not
#     ## be caught.
#     expect_warning(do.call(match.fun("check_heaping_whipple"),
#                            args = list(Value = 25:65, Age = 25:65)), NA)
#     expect_warning(do.call(check_heaping_whipple,
#                            args = list(Value = 25:65, Age = 25:65)), NA)
#     expect_warning(do.call(match.fun(check_heaping_whipple),
#                            args = list(Value = 25:65, Age = 25:65)), NA)
# 
# })
# 
# test_that("calling 'DemoTools::check_heaping_whipple' does not produce a warning", {
# 
#     expect_warning(DemoTools::check_heaping_whipple(Value = 25:65, Age = 25:65), NA)
#     expect_warning(DemoTools::check_heaping_whipple(Value = 25:65, Age = 25:65), NA)
# 
#     ## Invoking 'do.call' in these forms does not let 'match.call' see
#     ## the function name so any use of the old function name will not
#     ## be caught.
#     expect_warning(do.call(DemoTools::check_heaping_whipple,
#                            args = list(Value = 25:65, Age = 25:65)), NA)
#     expect_warning(do.call(match.fun(DemoTools::check_heaping_whipple),
#                            args = list(Value = 25:65, Age = 25:65)), NA)
# })
