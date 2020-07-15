#ifndef SIMPLINE_SIMPLINE_H
#define SIMPLINE_SIMPLINE_H

#include <Eigen/Dense>
#include <vector>
#include <map>

template<typename T>
struct simpline
{
	typedef typename Eigen::Matrix<T, 3, 1> Vector3;
	typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> VectorX;
	
	class ParametrizedSpline
	{
	public:
		ParametrizedSpline();
		
		ParametrizedSpline(const std::vector<T>& parameterValues, const std::vector<simpline<T>::Vector3>& points);
		
		simpline<T>::Vector3 getValue(const T& parameterValue) const;
		
		simpline<T>::Vector3 getGradient(const T& parameterValue) const;
		
		T getLength() const;
		
		T getLength(const T& startParameterValue, const T& endParameterValue) const;
	
	private:
		simpline<T>::Vector3 computeFiniteDifference(const size_t& startPointIndex, const size_t& endPointIndex);
		
		std::vector<T> parameterValues;
		std::vector<simpline<T>::Vector3> points;
		std::vector<simpline<T>::Vector3> firstDerivatives;
		std::vector<simpline<T>::Vector3> secondDerivatives;
		std::vector<simpline<T>::Vector3> thirdDerivatives;
		std::map<std::pair<size_t, size_t>, simpline<T>::Vector3> finiteDifferences;
	};
	
	class ConstantSpeedSpline
	{
	public:
		ConstantSpeedSpline();
		
		ConstantSpeedSpline(std::vector<simpline<T>::Vector3> points, T speed);
		
		simpline<T>::Vector3 getValue(const T& time) const;
		
		simpline<T>::Vector3 getGradient(const T& time) const;
		
		T getLength() const;
		
		T getLength(const T& startTime, const T& endTime) const;
		
		T getSpeed() const;
		
		T getDuration() const;
	
	private:
		T computeParameterValue(const T& time) const;
		
		std::map<T, T> timeParameterValues;
		ParametrizedSpline parametrizedSpline;
		T speed;
		T duration;
	};
};

#endif
