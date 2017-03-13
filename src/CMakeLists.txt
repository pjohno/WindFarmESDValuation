
SET(MY_LIBRARY_SRCS  
DelayedWindFarm.cpp 
delay_ops.cc            generic_SDE.cc       nodes.cc      param_maps.cc  simulator.cc       vary_params.cc
delay_ops_test.cc       generic_solver.cc    no_price.cc   param_vecs.cc  stencil.cc         vary_params_test.cc
diff_grid.cc            mtrand.cc            optimizer.cc  random_nd.cc   storage.cc
finite_diff_methods.cc  newton_iteration.cc  option.cc     sde.cc        

)

ADD_LIBRARY(DelayedWindFarmLIB SHARED
${MY_LIBRARY_SRCS} 
)

TARGET_LINK_LIBRARIES( DelayedWindFarmLIB boost_thread ${COMMONMATHLIB_LIBRARY})

ADD_LIBRARY(DelayedWindFarmLIB_static STATIC
${MY_LIBRARY_SRCS} 
)

TARGET_LINK_LIBRARIES( DelayedWindFarmLIB_static boost_thread ${COMMONMATHLIB_LIBRARY}_static)

