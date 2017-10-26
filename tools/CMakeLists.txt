include_directories(../src)

add_executable(index2ft index2ft.c)
target_link_libraries(index2ft fmftsaxs)

add_executable(single_saxs_fit single_saxs_fit.c)
target_link_libraries(single_saxs_fit fmftsaxs)

add_executable(puresaxs puresaxs.c)
target_link_libraries(puresaxs fmftsaxs)

add_executable(ft2euler ft2euler.c)
target_link_libraries(ft2euler fmftsaxs) 

add_executable(correlate correlate.c)
target_link_libraries(correlate fmftsaxs)

add_executable(score_ft_naive score_ft_naive.c)
target_link_libraries(score_ft_naive fmftsaxs)

add_executable(single_saxs single_saxs.c)
target_link_libraries(single_saxs fmftsaxs)
