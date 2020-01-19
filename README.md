# simpline
Simple constant-speed natural cubic spline interpolator for 3D points.

## Compilation
To compile this library and install it, in a terminal window, go to the root folder of the library and enter the following commands:
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
sudo make install
```

## Usage example
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
	simpline<float>::ConstantSpeedSpline spline(points, 0.5);
	
	float currentTime = 0.0;
	while(currentTime <= spline.getDuration())
	{
		std::cout << "Position (" << currentTime << "s):" << std::endl;
		std::cout << spline.getPosition(currentTime) << std::endl << std::endl;
		
		std::cout << "Gradient (" << currentTime << "s):" << std::endl;
		std::cout << spline.getGradient(currentTime) << std::endl << std::endl;
		
		currentTime += 1.0;
	}

	return 0;
}
```
