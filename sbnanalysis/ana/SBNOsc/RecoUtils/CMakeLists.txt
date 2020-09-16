add_library(sbnanalysis_SBNOscD_RecoUtils SHARED RecoUtils.cc)
include_directories($ENV{NUG4_INC})

target_link_libraries(sbnanalysis_SBNOscD_RecoUtils
  sbnanalysis_ProviderManager
  lardataobj_RawData
  lardata_RecoObjects
  lardataobj_RecoBase
)

install(TARGETS sbnanalysis_SBNOscD_RecoUtils DESTINATION lib)

