#include "Simpline.h"
#include "Constants.h"

template<typename T>
simpline<T>::ConstantSpeedSpline::ConstantSpeedSpline()
{
}

template<typename T>
simpline<T>::ConstantSpeedSpline::ConstantSpeedSpline(const std::vector<simpline<T>::Vector3>& points, const T& speed):
		timeParameterValues(), parametrizedSpline(), speed(speed), duration()
{
	if(points.size() < 2)
	{
		throw std::runtime_error("Point list must contain at least two items!");
	}
	
	if(speed <= 0)
	{
		throw std::runtime_error("Speed must be above 0!");
	}
	
	std::vector<T> parameterValues = { 0 };
	for(size_t i = 1; i < points.size(); i++)
	{
		parameterValues.push_back(parameterValues[i - 1] + (points[i] - points[i - 1]).norm());
	}
	parametrizedSpline = ParametrizedSpline(parameterValues, points);
	
	T splineLength = parametrizedSpline.getLength();
	duration = splineLength / speed;
	
	T time = 0;
	for(size_t i = 0; i < parameterValues.size() - 1; i++)
	{
		timeParameterValues[time] = parameterValues[i];
		time += (parametrizedSpline.getLength(parameterValues[i], parameterValues[i + 1]) / splineLength) * duration;
	}
	timeParameterValues[duration] = parameterValues[parameterValues.size() - 1];
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getValue(const T& time) const
{
	if(timeParameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get value from empty constant-speed spline. Use non-default constructor to provide points.");
	}
	
	if(time < 0.0 || time > duration)
	{
		throw std::runtime_error("Value requested at time=" + std::to_string(time) + ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}
	
	return parametrizedSpline.getValue(computeParameterValue(time));
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::computeParameterValue(const T& time) const
{
	if(timeParameterValues.find(time) != timeParameterValues.end())
	{
		return timeParameterValues.at(time);
	}
	
	const auto nextTimeParameterValue = timeParameterValues.upper_bound(time);
	const auto previousTimeParameterValue = std::prev(nextTimeParameterValue);
	const T wantedSplineLength = (speed * time) - parametrizedSpline.getLength(0, previousTimeParameterValue->second);
	
	if(wantedSplineLength <= 0)
	{
		return previousTimeParameterValue->second;
	}
	
	// bisection method
	T lowerBoundT = previousTimeParameterValue->second;
	T upperBoundT = nextTimeParameterValue->second;
	
	while(upperBoundT - lowerBoundT > epsilon)
	{
		const T middleT = (lowerBoundT + upperBoundT) / 2;
		const T lowerBoundSplineLength = parametrizedSpline.getLength(previousTimeParameterValue->second, lowerBoundT);
		const T middleSplineLength = parametrizedSpline.getLength(previousTimeParameterValue->second, middleT);
		
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
typename simpline<T>::Vector3 simpline<T>::ConstantSpeedSpline::getGradient(const T& time) const
{
	if(timeParameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get gradient from empty constant-speed spline. Use non-default constructor to provide points.");
	}
	
	if(time < 0.0 || time > duration)
	{
		throw std::runtime_error("Gradient requested at time=" + std::to_string(time) + ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}
	
	return parametrizedSpline.getGradient(computeParameterValue(time)).normalized() * speed;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getLength() const
{
	if(timeParameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get length of empty constant-speed spline. Use non-default constructor to provide points.");
	}
	
	return duration * speed;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getLength(const T& startTime, const T& endTime) const
{
	if(timeParameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get length of empty constant-speed spline. Use non-default constructor to provide points.");
	}
	
	if(startTime > endTime)
	{
		throw std::runtime_error("Length requested from time=" + std::to_string(startTime) + " to time=" + std::to_string(endTime) +
								 ". End time must be greater or equal to start time.");
	}
	
	if(startTime < 0.0 || endTime > duration)
	{
		throw std::runtime_error("Length requested from time=" + std::to_string(startTime) + " to time=" + std::to_string(endTime) +
								 ". Time must be between 0.0 and " + std::to_string(duration) + ".");
	}
	
	return (endTime - startTime) * speed;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getSpeed() const
{
	if(timeParameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get speed of empty constant-speed spline. Use non-default constructor to provide points.");
	}
	
	return speed;
}

template<typename T>
T simpline<T>::ConstantSpeedSpline::getDuration() const
{
	if(timeParameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get duration of empty constant-speed spline. Use non-default constructor to provide points.");
	}
	
	return duration;
}

template class simpline<float>::ConstantSpeedSpline;

template class simpline<double>::ConstantSpeedSpline;
