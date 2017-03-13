ADD_EXECUTABLE( DelayedWindFarm.out main.cpp )
 TARGET_LINK_LIBRARIES( DelayedWindFarm.out DelayedWindFarmLIB )

ADD_EXECUTABLE( Model03.out model03.cpp )
 TARGET_LINK_LIBRARIES( Model03.out DelayedWindFarmLIB )
