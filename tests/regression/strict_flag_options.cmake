# #45: suffixes on flags are errors, not silently selected solver settings.
foreach(option -gx -cx -rx -ix -vfoo -ex -ox -oix -apx -aPx -pnx -pNx -pjx -pJx)
    fastercap_test(trailing_flag${option} three_cubes 0 NO_MATRIX EXIT_CODE 64 ${option})
endforeach()
# Existing case variants remain accepted; no adaptive refinement is enabled here.
fastercap_test(preconditioner_flag_case three_cubes 0.01 SAME_AS -pn -pN)
fastercap_test(jacobi_flag_case three_cubes 0.01 SAME_AS -pj -pJ)
fastercap_test(auto_preconditioner_flag_case three_cubes 0.01 SAME_AS -ap -aP)
