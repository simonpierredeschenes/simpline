#include "Simpline.h"
#include "Constants.h"

template<typename T>
simpline<T>::ConstantSpeedSpline::ConstantSpeedSpline(std::vector<simpline<T>::Vector3> points, T speed):
		speed(speed)
{
	if(points.size() < 2)
	{
		throw std::runtime_error("Point list must contain at least two items!");
	}
	
	if(speed <= 0)
	{
		throw std::runtime_error("Speed must be above 0!");
	}
	
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
	
	const auto nextTimeDistance = timeDistances.upper_bound(time);
	const auto previousTimeDistance = std::prev(nextTimeDistance);
	const T wantedSplineLength = (speed * time) - parametrizedSpline.computePartLength(0, previousTimeDistance->second);
	
	if(wantedSplineLength <= 0)
	{
		return previousTimeDistance->second;
	}
	
	// bisection method
	T lowerBoundT = previousTimeDistance->second;
	T upperBoundT = nextTimeDistance->second;
	
	while(upperBoundT - lowerBoundT > epsilon)
	{
		const T middleT = (lowerBoundT + upperBoundT) / 2;
		const T lowerBoundSplineLength = parametrizedSpline.computePartLength(previousTimeDistance->second, lowerBoundT);
		const T middleSplineLength = parametrizedSpline.computePartLength(previousTimeDistance->second, middleT);
		
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
	if(time < 0.0 || time > duration)
	{
		throw std::runtime_error("Position requested at time=" + std::to_string(time) + ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}
	
	return parametrizedSpline.getPosition(computeDistance(time));
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getGradient(const T& time) const
{
	if(time < 0.0 || time > duration)
	{
		throw std::runtime_error("Gradient requested at time=" + std::to_string(time) + ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}
	
	return parametrizedSpline.getGradient(computeDistance(time)).normalized() * speed;
}

template class simpline<float>::ConstantSpeedSpline;

template class simpline<double>::ConstantSpeedSpline;
