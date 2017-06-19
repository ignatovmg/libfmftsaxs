include_directories(../src)

add_executable(ft2euler ft2euler.c)
target_link_libraries(ft2euler mfftsaxs) 

add_executable(fitmassha2_new fitmassha2_new.c)
target_link_libraries(fitmassha2_new mfftsaxs)

add_executable(mytest test.c)
target_link_libraries(mytest mfftsaxs)

add_executable(score_ft_silly score_ft_silly.c)
target_link_libraries(score_ft_silly mfftsaxs)
