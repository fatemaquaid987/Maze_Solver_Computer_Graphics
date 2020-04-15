# Steps to set up OpenGL for windows using Visual Studio.

The first step is to download Microsoft Visual Studio.  

Next, you will need to download certain OpenGL libraries which are listed are follows. Simply click on each of the library and download the zip files depending on your machine’s specifications(x64 or x86)  

1. [GLM](http://glm.g-truc.net/0.9.7/index.html) : Mathematics Library  
2. [GLFW](https://www.glfw.org/download.html) : creating windows  
3. [GLEW](http://glew.sourceforge.net/index.html) : linking other libraries  
4. [GLUT](http://freeglut.sourceforge.net/index.php#download) : system-specific tasks
  
In addition to above libraries you will need OpenGL library which comes pre-installed with Windows OS. (search the C drive for opengl32.lib)  

Next, make a folder where you will keep all OpenGL files. Make three subfolders as follows:  
1. headers : where you will keep all .h files (directly paste the contents of include folder)  
2. libraries : .lib files  
3. dlls : .dll files  

Unzip the downloaded library files and transfer the required files in corresponding sub-folders.  (Don’t forget to add the opengl.lib in the libraries folder)
For eg: .lib files in libraries folder and so on. 


The next step is to open Visual Studio. Follow the following steps in Visual Studio.
1. Go to file->New->project->empty project. Click Ok. A new project will be created.  
2. Right click on Source files in Solution explorer->add->new Item->C++ file.  
3. Right click on the project->Properties.  In Configuration -Select All configurations & in Platforms-select All Platforms.  
4. Select the VC++ directories->
Include Directories -> Edit ->new -> Browse -> (select the headers folder where you have kept all the header files.  
Library Directories -> Edit ->new -> Browse -> (select the libraries folder where you have kept all the library files.  
5. Select C/C++->General->Additional Include Directories -> Edit ->new -> Browse -> (select the headers folder where you have kept all the header files.  
6. Select Linker->General->Additional Library Directories -> Edit ->new -> Browse -> (select the libraries folder where you have kept all the library files.  
7. Select Linker->Input->Additional Dependencies-> Edit  
Add all the names of .lib files present in libraries folder.  
8. Right Click on the project-> Open Project in File Explorer -> copy & paste all the .dll files here.  

Now you are all ready to run your OpenGL program.  
	

Reference: http://www3.ntu.edu  
For additional help you can follow this [video](https://www.youtube.com/watch?v=tmTLcyAwBJo)   
