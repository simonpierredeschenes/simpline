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
		ParametrizedSpline(const std::vector<T>& ts, const std::vector<simpline<T>::Vector3>& points);
		simpline<T>::Vector3 getPosition(const T& t) const;
		simpline<T>::Vector3 getGradient(const T& t) const;
		T computePartLength(const T& t_start, const T& t_end) const;
	
	private:
		simpline<T>::Vector3 computeFiniteDifference(const size_t& startIndex, const size_t& endIndex);
		
		std::vector<T> ts;
		std::vector<simpline<T>::Vector3> points;
		std::vector<simpline<T>::Vector3> firstDerivatives;
		std::vector<simpline<T>::Vector3> secondDerivatives;
		std::vector<simpline<T>::Vector3> thirdDerivatives;
		std::map<std::pair<size_t, size_t>, simpline<T>::Vector3> finiteDifferences;
	};
	
	class ConstantSpeedSpline
	{
	public:
		ConstantSpeedSpline(std::vector<simpline<T>::Vector3> points, T speed);
		
		T getDuration() const;
		
		simpline<T>::Vector3 getPosition(const T& time) const;
		
		simpline<T>::Vector3 getGradient(const T& time) const;
	
	private:
		T computeDistance(const T& time) const;
		
		T duration;
		T speed;
		std::map<T, T> timeDistances;
		ParametrizedSpline parametrizedSpline;
	};
};


#endif
