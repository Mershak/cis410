# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# compile CXX with /usr/bin/c++
CXX_DEFINES = -DDIY_NO_THREADS -DFMT_SHARED -DH5_BUILT_AS_DYNAMIC_LIB -DVTK_HAS_OGGTHEORA_SUPPORT -Dkiss_fft_scalar=double -DvtkDomainsChemistry_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\" -DvtkFiltersCore_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\" -DvtkIOExportGL2PS_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\" -DvtkIOExport_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\" -DvtkRenderingContext2D_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\" -DvtkRenderingCore_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\" -DvtkRenderingOpenGL2_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\" -DvtkRenderingVolume_AUTOINIT_INCLUDE=\"/home/mert/Desktop/cis410/project2/CMakeFiles/vtkModuleAutoInit_bd64b765fc8236aa4bea0b628e677e8c.h\"

CXX_INCLUDES = -isystem /home/mert/Desktop/cis410/build/Wrapping/Tools -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Wrapping/Tools -isystem /home/mert/Desktop/cis410/build/Views/Infovis -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Views/Infovis -isystem /home/mert/Desktop/cis410/build/Common/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/Core -isystem /home/mert/Desktop/cis410/build/Utilities/KWIML -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Utilities/KWIML -isystem /home/mert/Desktop/cis410/build/Utilities/KWSys -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Utilities/KWSys -isystem /home/mert/Desktop/cis410/build/Common/DataModel -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/DataModel -isystem /home/mert/Desktop/cis410/build/Common/Math -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/Math -isystem /home/mert/Desktop/cis410/build/ThirdParty/kissfft/vtkkissfft -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/kissfft/vtkkissfft -isystem /home/mert/Desktop/cis410/build/ThirdParty/kissfft -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/kissfft -isystem /home/mert/Desktop/cis410/build/Common/Transforms -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/Transforms -isystem /home/mert/Desktop/cis410/build/Common/ExecutionModel -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/ExecutionModel -isystem /home/mert/Desktop/cis410/build/Interaction/Style -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Interaction/Style -isystem /home/mert/Desktop/cis410/build/Rendering/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/Core -isystem /home/mert/Desktop/cis410/build/Filters/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Core -isystem /home/mert/Desktop/cis410/build/Common/Misc -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/Misc -isystem /home/mert/Desktop/cis410/build/Rendering/Context2D -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/Context2D -isystem /home/mert/Desktop/cis410/build/Views/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Views/Core -isystem /home/mert/Desktop/cis410/build/Interaction/Widgets -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Interaction/Widgets -isystem /home/mert/Desktop/cis410/build/Filters/General -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/General -isystem /home/mert/Desktop/cis410/build/Filters/Sources -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Sources -isystem /home/mert/Desktop/cis410/build/Common/Color -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/Color -isystem /home/mert/Desktop/cis410/build/Views/Context2D -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Views/Context2D -isystem /home/mert/Desktop/cis410/build/ThirdParty/loguru/vtkloguru -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/loguru/vtkloguru -isystem /home/mert/Desktop/cis410/build/ThirdParty/loguru -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/loguru -isystem /home/mert/Desktop/cis410/build/Testing/Rendering -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Testing/Rendering -isystem /home/mert/Desktop/cis410/build/Testing/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Testing/Core -isystem /home/mert/Desktop/cis410/build/Rendering/VolumeOpenGL2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/VolumeOpenGL2 -isystem /home/mert/Desktop/cis410/build/Imaging/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Core -isystem /home/mert/Desktop/cis410/build/Imaging/Math -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Math -isystem /home/mert/Desktop/cis410/build/Rendering/OpenGL2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/OpenGL2 -isystem /home/mert/Desktop/cis410/build/Rendering/UI -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/UI -isystem /home/mert/Desktop/cis410/build/ThirdParty/glew/vtkglew -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/glew/vtkglew -isystem /home/mert/Desktop/cis410/build/ThirdParty/glew -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/glew -isystem /home/mert/Desktop/cis410/build/Rendering/Volume -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/Volume -isystem /home/mert/Desktop/cis410/build/Rendering/Label -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/Label -isystem /home/mert/Desktop/cis410/build/Rendering/FreeType -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/FreeType -isystem /home/mert/Desktop/cis410/build/ThirdParty/freetype/vtkfreetype -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/freetype/vtkfreetype -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/freetype/vtkfreetype/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/freetype/vtkfreetype/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/freetype -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/freetype -isystem /home/mert/Desktop/cis410/build/ThirdParty/zlib/vtkzlib -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/zlib/vtkzlib -isystem /home/mert/Desktop/cis410/build/ThirdParty/zlib -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/zlib -isystem /home/mert/Desktop/cis410/build/Utilities/octree -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Utilities/octree -isystem /home/mert/Desktop/cis410/build/Rendering/LOD -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/LOD -isystem /home/mert/Desktop/cis410/build/Rendering/Image -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/Image -isystem /home/mert/Desktop/cis410/build/Rendering/ContextOpenGL2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/ContextOpenGL2 -isystem /home/mert/Desktop/cis410/build/IO/VeraOut -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/VeraOut -isystem /home/mert/Desktop/cis410/build/ThirdParty/hdf5/vtkhdf5 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/hdf5/vtkhdf5 -isystem /home/mert/Desktop/cis410/build/ThirdParty/hdf5 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/hdf5 -isystem /home/mert/Desktop/cis410/build/ThirdParty/hdf5/vtkhdf5/src -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/hdf5/vtkhdf5/src -isystem /home/mert/Desktop/cis410/build/ThirdParty/hdf5/vtkhdf5/hl/src -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/hdf5/vtkhdf5/hl/src -isystem /home/mert/Desktop/cis410/build/IO/TecplotTable -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/TecplotTable -isystem /home/mert/Desktop/cis410/build/IO/SegY -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/SegY -isystem /home/mert/Desktop/cis410/build/IO/Image -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Image -isystem /home/mert/Desktop/cis410/build/IO/ParallelXML -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/ParallelXML -isystem /home/mert/Desktop/cis410/build/IO/XML -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/XML -isystem /home/mert/Desktop/cis410/build/IO/XMLParser -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/XMLParser -isystem /home/mert/Desktop/cis410/build/IO/PLY -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/PLY -isystem /home/mert/Desktop/cis410/build/IO/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Core -isystem /home/mert/Desktop/cis410/build/IO/OggTheora -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/OggTheora -isystem /home/mert/Desktop/cis410/build/IO/Movie -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Movie -isystem /home/mert/Desktop/cis410/build/ThirdParty/theora/vtktheora -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/theora/vtktheora -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/theora/vtktheora/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/theora -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/theora -isystem /home/mert/Desktop/cis410/build/ThirdParty/ogg/vtkogg -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/ogg/vtkogg -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/ogg/vtkogg/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/ogg/vtkogg/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/ogg -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/ogg -isystem /home/mert/Desktop/cis410/build/IO/NetCDF -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/NetCDF -isystem /home/mert/Desktop/cis410/build/ThirdParty/netcdf/vtknetcdf -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/netcdf/vtknetcdf -isystem /home/mert/Desktop/cis410/build/ThirdParty/netcdf -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/netcdf -isystem /home/mert/Desktop/cis410/build/IO/MotionFX -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/MotionFX -isystem /home/mert/Desktop/cis410/build/ThirdParty/pegtl -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/pegtl -isystem /home/mert/Desktop/cis410/build/IO/Parallel -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Parallel -isystem /home/mert/Desktop/cis410/build/IO/Geometry -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Geometry -isystem /home/mert/Desktop/cis410/build/IO/Legacy -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Legacy -isystem /home/mert/Desktop/cis410/build/ThirdParty/jsoncpp/vtkjsoncpp -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/jsoncpp/vtkjsoncpp -isystem /home/mert/Desktop/cis410/build/ThirdParty/jsoncpp/vtkjsoncpp/json -isystem /home/mert/Desktop/cis410/build/ThirdParty/jsoncpp -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/jsoncpp -isystem /home/mert/Desktop/cis410/build/IO/MINC -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/MINC -isystem /home/mert/Desktop/cis410/build/IO/LSDyna -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/LSDyna -isystem /home/mert/Desktop/cis410/build/IO/Infovis -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Infovis -isystem /home/mert/Desktop/cis410/build/ThirdParty/libxml2/vtklibxml2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libxml2/vtklibxml2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libxml2/vtklibxml2/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/libxml2/vtklibxml2/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/libxml2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libxml2 -isystem /home/mert/Desktop/cis410/build/IO/Import -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Import -isystem /home/mert/Desktop/cis410/build/IO/IOSS -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/IOSS -isystem /home/mert/Desktop/cis410/build/Parallel/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Parallel/Core -isystem /home/mert/Desktop/cis410/build/ThirdParty/ioss/vtkioss -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/ioss/vtkioss -isystem /home/mert/Desktop/cis410/build/ThirdParty/ioss -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/ioss -isystem /home/mert/Desktop/cis410/build/ThirdParty/exodusII/vtkexodusII -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/exodusII/vtkexodusII -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/exodusII/vtkexodusII/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/exodusII/vtkexodusII/include -isystem /home/mert/Desktop/cis410/build/ThirdParty/exodusII -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/exodusII -isystem /home/mert/Desktop/cis410/build/ThirdParty/cgns/vtkcgns/src -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/cgns/vtkcgns/src -isystem /home/mert/Desktop/cis410/build/ThirdParty/cgns -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/cgns -isystem /home/mert/Desktop/cis410/build/IO/HDF -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/HDF -isystem /home/mert/Desktop/cis410/build/IO/Video -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Video -isystem /home/mert/Desktop/cis410/build/IO/ExportPDF -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/ExportPDF -isystem /home/mert/Desktop/cis410/build/IO/Export -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Export -isystem /home/mert/Desktop/cis410/build/Rendering/VtkJS -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/VtkJS -isystem /home/mert/Desktop/cis410/build/Rendering/SceneGraph -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/SceneGraph -isystem /home/mert/Desktop/cis410/build/ThirdParty/libharu/vtklibharu/src -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libharu/vtklibharu/src -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libharu/vtklibharu/src/../include -isystem /home/mert/Desktop/cis410/build/ThirdParty/libharu/vtklibharu/src/../include -isystem /home/mert/Desktop/cis410/build/ThirdParty/libharu -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libharu -isystem /home/mert/Desktop/cis410/build/IO/ExportGL2PS -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/ExportGL2PS -isystem /home/mert/Desktop/cis410/build/Rendering/GL2PSOpenGL2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/GL2PSOpenGL2 -isystem /home/mert/Desktop/cis410/build/ThirdParty/gl2ps/vtkgl2ps -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/gl2ps/vtkgl2ps -isystem /home/mert/Desktop/cis410/build/ThirdParty/gl2ps -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/gl2ps -isystem /home/mert/Desktop/cis410/build/ThirdParty/png/vtkpng -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/png/vtkpng -isystem /home/mert/Desktop/cis410/build/ThirdParty/png -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/png -isystem /home/mert/Desktop/cis410/build/IO/Exodus -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Exodus -isystem /home/mert/Desktop/cis410/build/IO/EnSight -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/EnSight -isystem /home/mert/Desktop/cis410/build/IO/CityGML -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/CityGML -isystem /home/mert/Desktop/cis410/build/ThirdParty/pugixml/vtkpugixml -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/pugixml/vtkpugixml -isystem /home/mert/Desktop/cis410/build/ThirdParty/pugixml -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/pugixml -isystem /home/mert/Desktop/cis410/build/IO/Chemistry -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Chemistry -isystem /home/mert/Desktop/cis410/build/IO/CONVERGECFD -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/CONVERGECFD -isystem /home/mert/Desktop/cis410/build/IO/CGNS -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/CGNS -isystem /home/mert/Desktop/cis410/build/IO/Asynchronous -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/Asynchronous -isystem /home/mert/Desktop/cis410/build/IO/AMR -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/AMR -isystem /home/mert/Desktop/cis410/build/Interaction/Image -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Interaction/Image -isystem /home/mert/Desktop/cis410/build/Imaging/Stencil -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Stencil -isystem /home/mert/Desktop/cis410/build/Imaging/Statistics -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Statistics -isystem /home/mert/Desktop/cis410/build/Imaging/Morphological -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Morphological -isystem /home/mert/Desktop/cis410/build/Imaging/General -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/General -isystem /home/mert/Desktop/cis410/build/Imaging/Fourier -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Fourier -isystem /home/mert/Desktop/cis410/build/IO/SQL -isystem /home/mert/Desktop/cis410/VTK-9.1.0/IO/SQL -isystem /home/mert/Desktop/cis410/build/ThirdParty/sqlite/vtksqlite -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/sqlite/vtksqlite -isystem /home/mert/Desktop/cis410/build/ThirdParty/sqlite -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/sqlite -isystem /home/mert/Desktop/cis410/build/Geovis/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Geovis/Core -isystem /home/mert/Desktop/cis410/build/Infovis/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Infovis/Core -isystem /home/mert/Desktop/cis410/build/Imaging/Sources -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Sources -isystem /home/mert/Desktop/cis410/build/ThirdParty/libproj/vtklibproj/src -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libproj/vtklibproj/src -isystem /home/mert/Desktop/cis410/build/ThirdParty/libproj -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/libproj -isystem /home/mert/Desktop/cis410/build/Infovis/Layout -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Infovis/Layout -isystem /home/mert/Desktop/cis410/build/Rendering/Annotation -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Rendering/Annotation -isystem /home/mert/Desktop/cis410/build/Imaging/Hybrid -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Hybrid -isystem /home/mert/Desktop/cis410/build/Imaging/Color -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Imaging/Color -isystem /home/mert/Desktop/cis410/build/Filters/Topology -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Topology -isystem /home/mert/Desktop/cis410/build/Filters/Selection -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Selection -isystem /home/mert/Desktop/cis410/build/Filters/SMP -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/SMP -isystem /home/mert/Desktop/cis410/build/Filters/Programmable -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Programmable -isystem /home/mert/Desktop/cis410/build/Filters/Points -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Points -isystem /home/mert/Desktop/cis410/build/Filters/Modeling -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Modeling -isystem /home/mert/Desktop/cis410/build/Filters/Verdict -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Verdict -isystem /home/mert/Desktop/cis410/build/ThirdParty/verdict/vtkverdict -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/verdict/vtkverdict -isystem /home/mert/Desktop/cis410/build/ThirdParty/verdict -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/verdict -isystem /home/mert/Desktop/cis410/build/Filters/ParallelImaging -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/ParallelImaging -isystem /home/mert/Desktop/cis410/build/Filters/Imaging -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Imaging -isystem /home/mert/Desktop/cis410/build/Filters/Statistics -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Statistics -isystem /home/mert/Desktop/cis410/build/Filters/Parallel -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Parallel -isystem /home/mert/Desktop/cis410/build/Filters/Extraction -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Extraction -isystem /home/mert/Desktop/cis410/build/Filters/Geometry -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Geometry -isystem /home/mert/Desktop/cis410/build/Filters/Hybrid -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Hybrid -isystem /home/mert/Desktop/cis410/build/Filters/Texture -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Texture -isystem /home/mert/Desktop/cis410/build/Filters/HyperTree -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/HyperTree -isystem /home/mert/Desktop/cis410/build/Filters/Generic -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/Generic -isystem /home/mert/Desktop/cis410/build/Filters/FlowPaths -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/FlowPaths -isystem /home/mert/Desktop/cis410/build/Common/ComputationalGeometry -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/ComputationalGeometry -isystem /home/mert/Desktop/cis410/build/ThirdParty/eigen -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/eigen -isystem /home/mert/Desktop/cis410/build/Filters/AMR -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Filters/AMR -isystem /home/mert/Desktop/cis410/build/Domains/ChemistryOpenGL2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Domains/ChemistryOpenGL2 -isystem /home/mert/Desktop/cis410/build/Domains/Chemistry -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Domains/Chemistry -isystem /home/mert/Desktop/cis410/build/Charts/Core -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Charts/Core -isystem /home/mert/Desktop/cis410/build/Parallel/DIY -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Parallel/DIY -isystem /home/mert/Desktop/cis410/build/Common/System -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Common/System -isystem /home/mert/Desktop/cis410/build/ThirdParty/diy2 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/diy2 -isystem /home/mert/Desktop/cis410/build/ThirdParty/expat/vtkexpat -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/expat/vtkexpat -isystem /home/mert/Desktop/cis410/build/ThirdParty/expat -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/expat -isystem /home/mert/Desktop/cis410/build/ThirdParty/doubleconversion/vtkdoubleconversion -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/doubleconversion/vtkdoubleconversion -isystem /home/mert/Desktop/cis410/build/ThirdParty/doubleconversion -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/doubleconversion -isystem /home/mert/Desktop/cis410/build/ThirdParty/lz4/vtklz4 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/lz4/vtklz4 -isystem /home/mert/Desktop/cis410/build/ThirdParty/lz4 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/lz4 -isystem /home/mert/Desktop/cis410/build/ThirdParty/lzma/vtklzma -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/lzma/vtklzma -isystem /home/mert/Desktop/cis410/build/ThirdParty/lzma -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/lzma -isystem /home/mert/Desktop/cis410/build/ThirdParty/utf8 -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/utf8 -isystem /home/mert/Desktop/cis410/build/Utilities/DICOMParser -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Utilities/DICOMParser -isystem /home/mert/Desktop/cis410/build/ThirdParty/jpeg/vtkjpeg -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/jpeg/vtkjpeg -isystem /home/mert/Desktop/cis410/build/ThirdParty/jpeg -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/jpeg -isystem /home/mert/Desktop/cis410/build/Utilities/MetaIO/vtkmetaio -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Utilities/MetaIO/vtkmetaio -isystem /home/mert/Desktop/cis410/build/Utilities/MetaIO -isystem /home/mert/Desktop/cis410/VTK-9.1.0/Utilities/MetaIO -isystem /home/mert/Desktop/cis410/build/ThirdParty/tiff/vtktiff/libtiff -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/tiff/vtktiff/libtiff -isystem /home/mert/Desktop/cis410/build/ThirdParty/tiff -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/tiff -isystem /home/mert/Desktop/cis410/build/ThirdParty/fmt/vtkfmt -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/fmt/vtkfmt -isystem /home/mert/Desktop/cis410/build/ThirdParty/fmt -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/fmt -isystem /home/mert/Desktop/cis410/build/ThirdParty/exprtk -isystem /home/mert/Desktop/cis410/VTK-9.1.0/ThirdParty/exprtk

CXX_FLAGS = 

