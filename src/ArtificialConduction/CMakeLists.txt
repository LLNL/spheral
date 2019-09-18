set(ArtificialConduction_inst
    ArtificialConduction
    ArtificialConductionPolicy
    )

set(ArtificialConduction_sources )
instantiate(ArtificialConduction_inst ArtificialConduction_sources)

set(ArtificialConduction_headers
    ArtificialConduction.hh
    ArtificialConductionPolicy.hh
    ) 

blt_add_library(NAME	    Spheral_ArtificialConduction
		HEADERS     ${ArtificialConduction_headers}
		SOURCES     ${ArtificialConduction_sources}
		DEPENDS_ON  ${spheral_depends}
		)

install(TARGETS		    Spheral_ArtificialConduction
	EXPORT		    ${PROJECT_NAME}-targets
	DESTINATION	    lib
	INCLUDES DESTINATION include
	)

install(FILES	    ${ArtificialConduction_headers}
	DESTINATION include/ArtificialConduction
	)
