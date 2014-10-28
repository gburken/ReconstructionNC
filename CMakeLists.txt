PROJECT( reconstructionNC )
#https://github.com/sintefmath/sist/blob/master/CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
IF(COMMAND CMAKE_POLICY)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

INCLUDE ($ENV{BME_BASE_DIR}/Config/BMECMakeLists.txt)

LINK_DIRECTORIES ("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v6.0/lib/Win32") 

SET( putImagesIn3DFiles
	putImagesIn3D.cpp
	)
SET( putImagesIn3D_V2Files
	putImagesIn3D_V2.cpp
	)
SET( putImagesIn3D_PlusFiles
	putImagesIn3D_Plus.cpp
	)
# common stuff

INCLUDE_DIRECTORIES( $ENV{ALEX_COMMON_DIR} )

ADD_EXECUTABLE( putImagesIn3D ${putImagesIn3DFiles})
TARGET_LINK_LIBRARIES ( putImagesIn3D ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${OpenCV_LIBS} )

ADD_EXECUTABLE( putImagesIn3D_V2 ${putImagesIn3D_V2Files})
TARGET_LINK_LIBRARIES ( putImagesIn3D_V2 ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${OpenCV_LIBS} )

ADD_EXECUTABLE( putImagesIn3D_Plus ${putImagesIn3D_PlusFiles})
TARGET_LINK_LIBRARIES ( putImagesIn3D_Plus ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${OpenCV_LIBS} vtkPlusCommon )