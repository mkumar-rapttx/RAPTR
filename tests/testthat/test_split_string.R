context("String splitting")
library("RAPTR")

test_that("Changing delimiter",{
    expect_equal(split_string("abc.def.ghi",return_position = 1,sep_char = "."),"abc")
    expect_equal(split_string("abc_def_ghi",return_position = 1,sep_char = "_"),"abc")
    expect_equal(split_string("abc,def,ghi",return_position = 1,sep_char = ","),"abc")
})

test_that("Changing return position",{
    expect_equal(split_string("abc.def.ghi",return_position = 1,sep_char = "."),"abc")
    expect_equal(split_string("abc.def.ghi",return_position = 2,sep_char = "."),"def")
    expect_equal(split_string("abc.def.ghi",return_position = 3,sep_char = "."),"ghi")
})
