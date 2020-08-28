![](logo.png)

# simpline
simpline is a simple constant-speed natural cubic spline interpolation library for 3D points.

## Compilation
To compile this library and install it, in a terminal window, go to the root folder of the library and enter the following commands:
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
sudo make install
```

## Using simpline in CMake Projects
This library provides CMake support and can be used in one of your projects as in the following example:
```cmake
cmake_minimum_required(VERSION 2.8.3)
project(example_project)

find_package(simpline)

include_directories(${simpline_INCLUDE_DIRS})
add_executable(example example.cpp)
target_link_libraries(example ${simpline_LIBRARIES})
```

## Usage Example
Here is a code example of how the library can be used:
```c++
#include <simpline/Simpline.h>
#include <iostream>

const std::vector<Eigen::Vector3f> points = {
		Eigen::Vector3f(-1.0, -1.0, 0),
		Eigen::Vector3f(1.0, -1.0, 0.5),
		Eigen::Vector3f(1.0, 1.0, 1.0),
		Eigen::Vector3f(-1.0, 1.0, 1.5),
		Eigen::Vector3f(-1.0, -1.0, 2.0),
		Eigen::Vector3f(1.0, -1.0, 2.5),
		Eigen::Vector3f(1.0, 1.0, 3.0),
		Eigen::Vector3f(-1.0, 1.0, 3.5),
		Eigen::Vector3f(-1.0, -1.0, 4.0)
};

int main(int argc, char** argv)
{
	const float speed = 0.5;
	simpline<float>::ConstantSpeedSpline trajectory(points, speed);
	
	float currentTime = 0.0;
	while(currentTime <= trajectory.getDuration())
	{
		std::cout << "Position (" << currentTime << "s):" << std::endl;
		std::cout << trajectory.getValue(currentTime) << std::endl << std::endl;
		
		std::cout << "Velocity (" << currentTime << "s):" << std::endl;
		std::cout << trajectory.getGradient(currentTime) << std::endl << std::endl;
		
		currentTime += 1.0;
	}
	
	return 0;
}
```
