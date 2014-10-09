PROJECT( reconstructionNC )
#https://github.com/sintefmath/sist/blob/master/CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
IF(COMMAND CMAKE_POLICY)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)


# Microsoft specific items
IF( MSVC )
 # Kill the anoying MS VS warning about non-safe functions.
 # They hide real warnings.
 ADD_DEFINITIONS(
   /D_SCL_SECURE_NO_DEPRECATE
   /D_CRT_SECURE_NO_DEPRECATE
   /D_CRT_TIME_FUNCTIONS_NO_DEPRECATE
 )
ENDIF()

FIND_PACKAGE (ITK REQUIRED)
INCLUDE( ${USE_ITK_FILE} )

FIND_PACKAGE (VTK REQUIRED)
INCLUDE( ${USE_VTK_FILE} )

FIND_PACKAGE(OpenCV REQUIRED )

#LINK_LIBRARIES(${OpenCV_LIBS} )


#FIND_PACKAGE(CUDA REQUIRED)
 

LINK_DIRECTORIES ("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v6.0/lib/Win32") 


# Boost stuff

INCLUDE_DIRECTORIES( D:/VisualStudioBaseProjects/boost_1_56_0/bin32/include/boost-1_56)
LINK_DIRECTORIES( D:/VisualStudioBaseProjects/boost_1_56_0/bin32/Lib/${OutDir})

		

SET( putImagesIn3DFiles
	putImagesIn3D.cpp
	)
# common stuff

INCLUDE_DIRECTORIES( $ENV{ALEX_COMMON_DIR} )

ADD_EXECUTABLE( putImagesIn3D ${putImagesIn3DFiles})
TARGET_LINK_LIBRARIES ( putImagesIn3D ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${OpenCV_LIBS} )