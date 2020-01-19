#include "Simpline.h"
#include "Constants.h"

template<typename T>
simpline<T>::ConstantSpeedSpline::ConstantSpeedSpline(std::vector<simpline<T>::Vector3> points, T speed):
		speed(speed)
{
	std::vector<T> distances = { 0 };
	for(int i = 1; i < points.size(); i++)
	{
		distances.push_back(distances[i - 1] + (points[i] - points[i - 1]).norm());
	}
	parametrizedSpline = ParametrizedSpline(distances, points);
	
	T splineLength = parametrizedSpline.computePartLength(0, distances[distances.size() - 1]);
	duration = splineLength / speed;
	
	T time = 0;
	for(size_t i = 0; i < distances.size() - 1; i++)
	{
		timeDistances[time] = distances[i];
		time += (parametrizedSpline.computePartLength(distances[i], distances[i + 1]) / splineLength) * duration;
	}
	timeDistances[duration] = distances[distances.size() - 1];
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getDuration() const
{
	return duration;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::computeDistance(const T& time) const
{
	if(timeDistances.find(time) != timeDistances.end())
	{
		return timeDistances.at(time);
	}
	
	const T wantedSplineLength = speed * time;
	const auto nextTimeDistance = timeDistances.upper_bound(time);
	const auto previousTimeDistance = std::prev(nextTimeDistance);
	
	// bisection method
	T lowerBoundT = previousTimeDistance->second;
	T upperBoundT = nextTimeDistance->second;
	
	while(upperBoundT - lowerBoundT > epsilon)
	{
		const T middleT = (lowerBoundT + upperBoundT) / 2;
		const T lowerBoundSplineLength = parametrizedSpline.computePartLength(0, lowerBoundT); // only compute partial spline length?
		const T middleSplineLength = parametrizedSpline.computePartLength(0, middleT);
		
		if((lowerBoundSplineLength - wantedSplineLength) * (middleSplineLength - wantedSplineLength) < 0)
		{
			upperBoundT = middleT;
		}
		else
		{
			lowerBoundT = middleT;
		}
	}
	
	return (lowerBoundT + upperBoundT) / 2;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getPosition(const T& time) const
{
	if(time > duration)
	{
		throw std::runtime_error("Requested time " + std::to_string(time) + " is too high! Max time is " + std::to_string(duration) + ".");
	}
	
	return parametrizedSpline.getPosition(computeDistance(time));
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getGradient(const T& time) const
{
	if(time > duration)
	{
		throw std::runtime_error("Requested time " + std::to_string(time) + " is too high! Max time is " + std::to_string(duration) + ".");
	}
	
	return parametrizedSpline.getGradient(computeDistance(time));
}

template class simpline<float>::ConstantSpeedSpline;
template class simpline<double>::ConstantSpeedSpline;
