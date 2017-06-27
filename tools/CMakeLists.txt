include_directories(../src)

add_executable(ft2euler ft2euler.c)
target_link_libraries(ft2euler fmftsaxs) 

add_executable(correlate correlate.c)
target_link_libraries(correlate fmftsaxs)

add_executable(score_ft_naive score_ft_naive.c)
target_link_libraries(score_ft_naive fmftsaxs)

add_executable(single_saxs single_saxs.c)
target_link_libraries(single_saxs fmftsaxs)
