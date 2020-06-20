#include "Simpline.h"
#include "Constants.h"

template<typename T>
simpline<T>::ParametrizedSpline::ParametrizedSpline()
{
}

template<typename T>
simpline<T>::ParametrizedSpline::ParametrizedSpline(const std::vector<T>& ts, const std::vector<simpline<T>::Vector3>& points):
		ts(ts), points(points), firstDerivatives(points.size() - 1), secondDerivatives(points.size()), thirdDerivatives(points.size() - 1), finiteDifferences()
{
	if(points.size() < 2)
	{
		throw std::runtime_error("Point list must contain at least two items!");
	}
	
	if(ts.size() != points.size())
	{
		throw std::runtime_error("Number of ts must be equal to number of points!");
	}
	
	std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> A(3);
	for(size_t i = 0; i < 3; i++)
	{
		A[i] = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(points.size(), points.size());
		for(int j = 0; j < points.size() - 2; j++)
		{
			for(int k = 0; k < points.size(); k++)
			{
				T coefficient = 0.0;
				if((j + 1) == k)
				{
					coefficient = 2.0;
				}
				else if((j + 1) - k == -1)
				{
					coefficient = (ts[j + 2] - ts[j + 1]) / (ts[j + 2] - ts[j]);
				}
				else if((j + 1) - k == 1)
				{
					coefficient = (ts[j + 1] - ts[j]) / (ts[j + 2] - ts[j]);
				}
				
				A[i](j + 1, k) = coefficient;
			}
		}
		
		// second derivative at first point is known
		A[i](0, 0) = 1;
		for(size_t j = 1; j < points.size(); j++)
		{
			A[i](0, j) = 0;
		}
		
		// second derivative at last point is known
		A[i](points.size() - 1, points.size() - 1) = 1;
		for(size_t j = 0; j < points.size() - 1; j++)
		{
			A[i](points.size() - 1, j) = 0;
		}
	}
	
	std::vector<simpline<T>::VectorX> b(3);
	for(size_t i = 0; i < 3; i++)
	{
		b[i] = simpline<T>::VectorX(points.size());
		for(size_t j = 0; j < points.size() - 2; j++)
		{
			b[i](j + 1) = 6 * computeFiniteDifference(j, j + 2)[i];
		}
		
		// second derivative at first point is zero
		b[i](0) = 0;
		// second derivative at last point is zero
		b[i](points.size() - 1) = 0;
	}
	
	for(size_t i = 0; i < 3; i++)
	{
		simpline<T>::VectorX x(A[i].colPivHouseholderQr().solve(b[i]));
		for(size_t j = 0; j < points.size(); j++)
		{
			secondDerivatives[j](i) = x(j);
		}
	}
	
	for(size_t j = 0; j < points.size() - 1; j++)
	{
		firstDerivatives[j] = computeFiniteDifference(j, j + 1) -
							  ((ts[j + 1] - ts[j]) * secondDerivatives[j] / 3.0) -
							  ((ts[j + 1] - ts[j]) * secondDerivatives[j + 1] / 6.0);
		thirdDerivatives[j] = (secondDerivatives[j + 1] - secondDerivatives[j]) / (ts[j + 1] - ts[j]);
	}
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::computeFiniteDifference(const size_t& startIndex, const size_t& endIndex)
{
	std::pair<size_t, size_t> indexPair(startIndex, endIndex);
	if(finiteDifferences.find(indexPair) != finiteDifferences.end())
	{
		return finiteDifferences[indexPair];
	}
	
	simpline<T>::Vector3 finiteDifference;
	if(startIndex == endIndex)
	{
		finiteDifference = points[startIndex];
	}
	else
	{
		finiteDifference =
				(computeFiniteDifference(startIndex + 1, endIndex) - computeFiniteDifference(startIndex, endIndex - 1)) / (ts[endIndex] - ts[startIndex]);
	}
	
	finiteDifferences[indexPair] = finiteDifference;
	
	return finiteDifference;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getPosition(const T& t) const
{
	if(t < ts[0] || t > ts[ts.size() - 1])
	{
		throw std::runtime_error("Position requested at t=" + std::to_string(t) + ". t must be between " + std::to_string(ts[0]) + " and " +
								 std::to_string(ts[ts.size() - 1]) + ".");
	}
	
	int previousIndex = 0;
	while(previousIndex < ts.size() - 1 && ts[previousIndex + 1] <= t)
	{
		previousIndex++;
	}
	
	return points[previousIndex] +
		   firstDerivatives[previousIndex] * (t - ts[previousIndex]) +
		   (secondDerivatives[previousIndex] / 2.0) * std::pow((t - ts[previousIndex]), 2) +
		   (thirdDerivatives[previousIndex] / 6.0) * std::pow((t - ts[previousIndex]), 3);
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getGradient(const T& t) const
{
	if(t < ts[0] || t > ts[ts.size() - 1])
	{
		throw std::runtime_error("Gradient requested at t=" + std::to_string(t) + ". t must be between " + std::to_string(ts[0]) + " and " +
								 std::to_string(ts[ts.size() - 1]) + ".");
	}
	
	int previousIndex = 0;
	while(previousIndex < ts.size() - 1 && ts[previousIndex + 1] <= t)
	{
		previousIndex++;
	}
	
	return firstDerivatives[previousIndex] +
		   secondDerivatives[previousIndex] * (t - ts[previousIndex]) +
		   (thirdDerivatives[previousIndex] / 2.0) * std::pow((t - ts[previousIndex]), 2);
}

template<typename T>
T simpline<T>::ParametrizedSpline::computePartLength(const T& t_start, const T& t_end) const
{
	if(t_start > t_end)
	{
		throw std::runtime_error("Part length requested from t=" + std::to_string(t_start) + " to t=" + std::to_string(t_end) +
								 ". End time must be greater or equal to start time.");
	}
	
	if(t_start < ts[0] || t_end > ts[ts.size() - 1])
	{
		throw std::runtime_error("Part length requested from t=" + std::to_string(t_start) + " to t=" + std::to_string(t_end) + ". t must be between " +
								 std::to_string(ts[0]) + " and " + std::to_string(ts[ts.size() - 1]) + ".");
	}
	
	int firstPointIndex = 0;
	while(firstPointIndex < ts.size() - 1 && ts[firstPointIndex] < t_start)
	{
		firstPointIndex++;
	}
	
	int lastPointIndex = 0;
	while(lastPointIndex < ts.size() - 1 && ts[lastPointIndex + 1] <= t_end)
	{
		lastPointIndex++;
	}
	
	// Gaussian quadrature spline length computation, taken from https://medium.com/@all2one/how-to-compute-the-length-of-a-spline-e44f5f04c40
	T length = 0.0;
	
	// special case for when t_start and t_end are in the same spline segment
	if(firstPointIndex > lastPointIndex)
	{
		T intervalLength = t_end - t_start;
		for(size_t i = 0; i < gaussianQuadratureAbcissa.size(); i++)
		{
			const T t = t_start + (((gaussianQuadratureAbcissa[i] + 1.0) / 2.0) * intervalLength); // Change of interval from [-1, 1]
			length += (intervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[i]; // Same for (intervalLength / 2.0)
		}
		return length;
	}
	
	// start partial spline segment computation
	T startIntervalLength = ts[firstPointIndex] - t_start;
	for(size_t i = 0; i < gaussianQuadratureAbcissa.size(); i++)
	{
		const T t = t_start + (((gaussianQuadratureAbcissa[i] + 1.0) / 2.0) * startIntervalLength); // Change of interval from [-1, 1]
		length += (startIntervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[i]; // Same for (startIntervalLength / 2.0)
	}
	// full spline segments computation
	for(size_t i = firstPointIndex; i < lastPointIndex; i++)
	{
		T intervalLength = ts[i + 1] - ts[i];
		for(size_t j = 0; j < gaussianQuadratureAbcissa.size(); j++)
		{
			const T t = ts[i] + (((gaussianQuadratureAbcissa[j] + 1.0) / 2.0) * intervalLength); // Change of interval from [-1, 1]
			length += (intervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[j]; // Same for (intervalLength / 2.0)
		}
	}
	// end partial spline segment computation
	T endIntervalLength = t_end - ts[lastPointIndex];
	for(size_t i = 0; i < gaussianQuadratureAbcissa.size(); i++)
	{
		const T t = ts[lastPointIndex] + (((gaussianQuadratureAbcissa[i] + 1.0) / 2.0) * endIntervalLength); // Change of interval from [-1, 1]
		length += (endIntervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[i]; // Same for (endIntervalLength / 2.0)
	}
	
	return length;
}

template class simpline<float>::ParametrizedSpline;

template class simpline<double>::ParametrizedSpline;
