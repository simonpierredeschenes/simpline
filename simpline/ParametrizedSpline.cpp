#include "Simpline.h"
#include "Constants.h"
#include <numeric>

// taken from https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template<typename T>
std::vector<size_t> sortIndices(const std::vector<T>& unsortedVector)
{
	std::vector<size_t> sortedIndices(unsortedVector.size());
	std::iota(sortedIndices.begin(), sortedIndices.end(), 0);
	std::sort(sortedIndices.begin(), sortedIndices.end(), [&unsortedVector](size_t i1, size_t i2){ return unsortedVector[i1] < unsortedVector[i2]; });
	return sortedIndices;
}

template<typename T>
simpline<T>::ParametrizedSpline::ParametrizedSpline()
{
}

template<typename T>
simpline<T>::ParametrizedSpline::ParametrizedSpline(const std::vector<T>& parameterValues, const std::vector<simpline<T>::Vector3>& points):
		parameterValues(parameterValues.size()), points(points.size()), firstDerivatives(points.size() - 1), secondDerivatives(points.size()),
		thirdDerivatives(points.size() - 1), finiteDifferences()
{
	if(points.size() < 2)
	{
		throw std::runtime_error("Point list must contain at least two items!");
	}
	
	if(parameterValues.size() != points.size())
	{
		throw std::runtime_error("Number of parameter values must be equal to number of points!");
	}
	
	std::vector<size_t> sortedIndices = sortIndices(parameterValues);
	for(size_t i = 0; i < sortedIndices.size(); i++)
	{
		this->parameterValues[i] = parameterValues[sortedIndices[i]];
		this->points[i] = points[sortedIndices[i]];
		
		if(i > 0 && this->parameterValues[i - 1] == this->parameterValues[i])
		{
			throw std::runtime_error("Multiple points cannot have the same parameter value.");
		}
	}
	
	std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> A(3);
	for(size_t i = 0; i < 3; i++)
	{
		A[i] = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(this->points.size(), this->points.size());
		for(size_t j = 0; j < this->points.size() - 2; j++)
		{
			for(size_t k = 0; k < this->points.size(); k++)
			{
				T coefficient = 0.0;
				if((j + 1) == k)
				{
					coefficient = 2.0;
				}
				else if((j + 1) - k == -1)
				{
					coefficient = (this->parameterValues[j + 2] - this->parameterValues[j + 1]) / (this->parameterValues[j + 2] - this->parameterValues[j]);
				}
				else if((j + 1) - k == 1)
				{
					coefficient = (this->parameterValues[j + 1] - this->parameterValues[j]) / (this->parameterValues[j + 2] - this->parameterValues[j]);
				}
				
				A[i](j + 1, k) = coefficient;
			}
		}
		
		// second derivative at first point is known
		A[i](0, 0) = 1;
		for(size_t j = 1; j < this->points.size(); j++)
		{
			A[i](0, j) = 0;
		}
		
		// second derivative at last point is known
		A[i](this->points.size() - 1, this->points.size() - 1) = 1;
		for(size_t j = 0; j < this->points.size() - 1; j++)
		{
			A[i](this->points.size() - 1, j) = 0;
		}
	}
	
	std::vector<simpline<T>::VectorX> b(3);
	for(size_t i = 0; i < 3; i++)
	{
		b[i] = simpline<T>::VectorX(this->points.size());
		for(size_t j = 0; j < this->points.size() - 2; j++)
		{
			b[i](j + 1) = 6 * computeFiniteDifference(j, j + 2)[i];
		}
		
		// second derivative at first point is zero
		b[i](0) = 0;
		// second derivative at last point is zero
		b[i](this->points.size() - 1) = 0;
	}
	
	for(size_t i = 0; i < 3; i++)
	{
		simpline<T>::VectorX x(A[i].colPivHouseholderQr().solve(b[i]));
		for(size_t j = 0; j < this->points.size(); j++)
		{
			secondDerivatives[j](i) = x(j);
		}
	}
	
	for(size_t j = 0; j < this->points.size() - 1; j++)
	{
		firstDerivatives[j] = computeFiniteDifference(j, j + 1) -
							  ((this->parameterValues[j + 1] - this->parameterValues[j]) * secondDerivatives[j] / 3.0) -
							  ((this->parameterValues[j + 1] - this->parameterValues[j]) * secondDerivatives[j + 1] / 6.0);
		thirdDerivatives[j] = (secondDerivatives[j + 1] - secondDerivatives[j]) / (this->parameterValues[j + 1] - this->parameterValues[j]);
	}
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::computeFiniteDifference(const size_t& startPointIndex, const size_t& endPointIndex)
{
	std::pair<size_t, size_t> indexPair(startPointIndex, endPointIndex);
	if(finiteDifferences.find(indexPair) != finiteDifferences.end())
	{
		return finiteDifferences[indexPair];
	}
	
	simpline<T>::Vector3 finiteDifference;
	if(startPointIndex == endPointIndex)
	{
		finiteDifference = points[startPointIndex];
	}
	else
	{
		finiteDifference =
				(computeFiniteDifference(startPointIndex + 1, endPointIndex) - computeFiniteDifference(startPointIndex, endPointIndex - 1)) /
				(parameterValues[endPointIndex] - parameterValues[startPointIndex]);
	}
	
	finiteDifferences[indexPair] = finiteDifference;
	
	return finiteDifference;
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getValue(const T& parameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get value from empty parametrized spline. Use non-default constructor to provide points.");
	}
	
	if(parameterValue < parameterValues[0] || parameterValue > parameterValues[parameterValues.size() - 1])
	{
		throw std::runtime_error("Value requested at parameterValue=" + std::to_string(parameterValue) + ". Parameter Value must be between " +
								 std::to_string(parameterValues[0]) + " and " + std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
	}
	
	int previousPointIndex = 0;
	while(previousPointIndex < parameterValues.size() - 1 && parameterValues[previousPointIndex + 1] <= parameterValue)
	{
		previousPointIndex++;
	}
	
	return points[previousPointIndex] +
		   firstDerivatives[previousPointIndex] * (parameterValue - parameterValues[previousPointIndex]) +
		   (secondDerivatives[previousPointIndex] / 2.0) * std::pow((parameterValue - parameterValues[previousPointIndex]), 2) +
		   (thirdDerivatives[previousPointIndex] / 6.0) * std::pow((parameterValue - parameterValues[previousPointIndex]), 3);
}

template<typename T>
typename simpline<T>::Vector3 simpline<T>::ParametrizedSpline::getGradient(const T& parameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get gradient from empty parametrized spline. Use non-default constructor to provide points.");
	}
	
	if(parameterValue < parameterValues[0] || parameterValue > parameterValues[parameterValues.size() - 1])
	{
		throw std::runtime_error("Gradient requested at parameterValue=" + std::to_string(parameterValue) + ". Parameter Value must be between " +
								 std::to_string(parameterValues[0]) + " and " + std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
	}
	
	int previousPointIndex = 0;
	while(previousPointIndex < parameterValues.size() - 1 && parameterValues[previousPointIndex + 1] <= parameterValue)
	{
		previousPointIndex++;
	}
	
	return firstDerivatives[previousPointIndex] +
		   secondDerivatives[previousPointIndex] * (parameterValue - parameterValues[previousPointIndex]) +
		   (thirdDerivatives[previousPointIndex] / 2.0) * std::pow((parameterValue - parameterValues[previousPointIndex]), 2);
}

template<typename T>
T simpline<T>::ParametrizedSpline::getLength() const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get length of empty parametrized spline. Use non-default constructor to provide points.");
	}
	
	// Gaussian quadrature spline length computation, taken from https://medium.com/@all2one/how-to-compute-the-length-of-a-spline-e44f5f04c40
	T length = 0.0;
	for(size_t i = 0; i < parameterValues.size() - 1; i++)
	{
		T intervalLength = parameterValues[i + 1] - parameterValues[i];
		for(size_t j = 0; j < gaussianQuadratureAbcissa.size(); j++)
		{
			const T t = parameterValues[i] + (((gaussianQuadratureAbcissa[j] + 1.0) / 2.0) * intervalLength); // Change of interval from [-1, 1]
			length += (intervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[j]; // Same for (intervalLength / 2.0)
		}
	}
	
	return length;
}

template<typename T>
T simpline<T>::ParametrizedSpline::getLength(const T& startParameterValue, const T& endParameterValue) const
{
	if(parameterValues.size() == 0)
	{
		throw std::runtime_error("Cannot get length of empty parametrized spline. Use non-default constructor to provide points.");
	}
	
	if(startParameterValue > endParameterValue)
	{
		throw std::runtime_error("Length requested from parameterValue=" + std::to_string(startParameterValue) + " to parameterValue=" +
								 std::to_string(endParameterValue) + ". End parameter value must be greater or equal to start parameter value.");
	}
	
	if(startParameterValue < parameterValues[0] || endParameterValue > parameterValues[parameterValues.size() - 1])
	{
		throw std::runtime_error("Length requested from parameterValue=" + std::to_string(startParameterValue) + " to parameterValue=" +
								 std::to_string(endParameterValue) + ". Parameter value must be between " + std::to_string(parameterValues[0]) + " and " +
								 std::to_string(parameterValues[parameterValues.size() - 1]) + ".");
	}
	
	// index computation of first point within parameter values
	int firstPointIndex = 0;
	while(firstPointIndex < parameterValues.size() - 1 && parameterValues[firstPointIndex] < startParameterValue)
	{
		firstPointIndex++;
	}
	
	// index computation of last point within parameter values
	int lastPointIndex = 0;
	while(lastPointIndex < parameterValues.size() - 1 && parameterValues[lastPointIndex + 1] <= endParameterValue)
	{
		lastPointIndex++;
	}
	
	// Gaussian quadrature spline length computation, taken from https://medium.com/@all2one/how-to-compute-the-length-of-a-spline-e44f5f04c40
	T length = 0.0;
	
	// special case: length computation for when startParameterValue and endParameterValue are in the same spline segment
	if(firstPointIndex > lastPointIndex)
	{
		T intervalLength = endParameterValue - startParameterValue;
		for(size_t i = 0; i < gaussianQuadratureAbcissa.size(); i++)
		{
			const T t = startParameterValue + (((gaussianQuadratureAbcissa[i] + 1.0) / 2.0) * intervalLength); // Change of interval from [-1, 1]
			length += (intervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[i]; // Same for (intervalLength / 2.0)
		}
		return length;
	}
	
	// length computation of first spline segment (partial segment)
	T startIntervalLength = parameterValues[firstPointIndex] - startParameterValue;
	for(size_t i = 0; i < gaussianQuadratureAbcissa.size(); i++)
	{
		const T t = startParameterValue + (((gaussianQuadratureAbcissa[i] + 1.0) / 2.0) * startIntervalLength); // Change of interval from [-1, 1]
		length += (startIntervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[i]; // Same for (startIntervalLength / 2.0)
	}
	// length computation of full spline segments
	for(size_t i = firstPointIndex; i < lastPointIndex; i++)
	{
		T intervalLength = parameterValues[i + 1] - parameterValues[i];
		for(size_t j = 0; j < gaussianQuadratureAbcissa.size(); j++)
		{
			const T t = parameterValues[i] + (((gaussianQuadratureAbcissa[j] + 1.0) / 2.0) * intervalLength); // Change of interval from [-1, 1]
			length += (intervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[j]; // Same for (intervalLength / 2.0)
		}
	}
	// length computation of last spline segment (partial segment)
	T endIntervalLength = endParameterValue - parameterValues[lastPointIndex];
	for(size_t i = 0; i < gaussianQuadratureAbcissa.size(); i++)
	{
		const T t = parameterValues[lastPointIndex] + (((gaussianQuadratureAbcissa[i] + 1.0) / 2.0) * endIntervalLength); // Change of interval from [-1, 1]
		length += (endIntervalLength / 2.0) * getGradient(t).norm() * gaussianQuadratureWeights[i]; // Same for (endIntervalLength / 2.0)
	}
	
	return length;
}

template class simpline<float>::ParametrizedSpline;

template class simpline<double>::ParametrizedSpline;
