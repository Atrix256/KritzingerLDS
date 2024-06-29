#include <stdio.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <direct.h>

// 0 means randomize
#define RANDOM_SEED() 1337

static const float c_goldenRatio = 1.61803398874989484820f;
static const float c_goldenRatioConjugate = 0.61803398874989484820f;
static const float c_pi = 3.14159265359f;

float Lerp(float a, float b, float t)
{
	return a * (1.0f - t) + b * t;
}

std::mt19937 GetRNG()
{
	std::random_device rd;
	unsigned int seed = (RANDOM_SEED() == 0) ? rd() : RANDOM_SEED();
	std::mt19937 rng(seed);
	return rng;
}

void KritzingerLDS_AddPoint(std::vector<float>& p, std::vector<float>& sortedp)
{
	// Kritzinger Low Discrepancy Sequence
	// A sequence that outperforms the golden ratio LDS.
	// https://arxiv.org/abs/2406.18132
	//
	// If we have N points in p, we want to find a y value that minimizes the function F(y,p).
	// F(y,p) = (N+1)*y^2 - y - 2*Sum(max(x,y))
	// the Sum is for all x in p.
	// 
	// If we have N points sorted x_1 < ... < x_n, also with x_0 = 0 and x_(n+1) = 1
	// then we have N=1 regions between points.
	// We can solve the function at each region and find whichever is the smallest.
	//
	// To find the minimum of the quadratic, we differentiate, and find where that equals 0.
	//
	// F(y,p) = Ay^2 + By + C
	// F'(y,p) = 2Ay + B
	// 
	// To find the minimum point of F, we find where F' is 0.
	// 2Ay+B = 0 : 2Ay = -B : y = -B/(2A)
	//
	// If y isn't in range of our section, we ignore it.
	// We could clamp it instead of ignoring, but we already have the end points in our point list.
	// If it is in range, we evaluate F at that point and keep the
	// y value that is smallest for all regions.
	// That is our new point

	const size_t N = p.size();

	const float A = float(N + 1);

	float bestY = 0.0f;
	float bestYScore = FLT_MAX;

	for (size_t segmentIndex = 0; segmentIndex <= N; ++segmentIndex)
	{
		// calculate the sum term
		float B = -1.0f;
		float C = 0.0;
		for (size_t xIndex = 0; xIndex < N; ++xIndex)
		{
			if (xIndex < segmentIndex)
				B -= 2.0f;
			else
				C -= 2.0f * sortedp[xIndex];
		}

		// get the edges of our region
		float lastx = (segmentIndex > 0) ? sortedp[segmentIndex - 1] : 0.0f;
		float x = (segmentIndex < N) ? sortedp[segmentIndex] : 1.0f;

		// calculate the minimum y and it's score
		// Verify y is in range
		float y = -B / (2.0f * A);
		if (y <= lastx || y >= x)
			continue;
		float yScore = A * y * y + B * y + C;
		
		// Keep the best scoring y
		if (yScore < bestYScore)
		{
			bestYScore = yScore;
			bestY = y;
		}
	}

	// Add the point to the list
	p.push_back(bestY);
	
	// Add the point to the sorted list
	sortedp.push_back(bestY);
	std::sort(sortedp.begin(), sortedp.end());
}

void KritzingerLDS(std::vector<float>& p, size_t count)
{
	std::vector<float> sortedp;
	p.resize(0);
	p.reserve(count);
	for (size_t i = 0; i < count; ++i)
		KritzingerLDS_AddPoint(p, sortedp);
}

void RegularPoints(std::vector<float>& p, size_t count)
{
	p.resize(count);
	for (size_t i = 0; i < count; ++i)
		p[i] = (float(i) + 0.5f) / float(count);
}

void RegularPointsTouchWalls(std::vector<float>& p, size_t count)
{
	p.resize(count);

	if (count == 1)
	{
		p[0] = 0.5f;
		return;
	}

	for (size_t i = 0; i < count; ++i)
		p[i] = float(i) / float(count-1);
}

void GoldenRatioLDS(std::vector<float>& p, size_t count)
{
	float value = 0.0f;
	p.resize(count);
	for (size_t i = 0; i < count; ++i)
	{
		p[i] = value;
		value = std::fmodf(value + c_goldenRatioConjugate, 1.0f);
	}
}

void Stratified(std::vector<float>& p, size_t count)
{
	std::mt19937 rng = GetRNG();
	std::uniform_real_distribution<float> dist(0.0f, 1.0f);

	p.resize(count);
	for (size_t i = 0; i < count; ++i)
		p[i] = (float(i) + dist(rng)) / float(count);
}

float Function_Triangle(float x)
{
	// Integrating this from 0 to 1 gives 0.5
	return x;
}

float Function_Step(float x)
{
	// Integrating this from 0 to 1 gives 0.6
	return (x < 0.6f) ? 1.0f : 0.0f;
}

float Function_Sine(float x)
{
	// Integrating this from 0 to 1 gives 2 / pi.
	return std::sin(x * c_pi);
}

struct Result
{
	const char* name = nullptr;
	std::vector<float> values;
};

template <typename TIntegrationFn, typename TPointFn>
Result IntegrationTestPoints_Sequence(const TIntegrationFn& IntegrationFn, float actualValue, const char* name, const TPointFn& PointFn, size_t pointCount)
{
	printf("  %s\n", name);

	// prepare our results
	Result ret;
	ret.name = name;
	ret.values.resize(pointCount);

	// Generate points
	std::vector<float> points;
	PointFn(points, pointCount);

	// Integerate
	float avgY = 0.0f;
	for (size_t index = 0; index < pointCount; ++index)
	{
		float y = IntegrationFn(points[index]);
		avgY = Lerp(avgY, y, 1.0f / float(index + 1));

		float error = avgY - actualValue;

		ret.values[index] = std::abs(error);
	}

	// return result
	return ret;
}

template <typename TIntegrationFn, typename TPointFn>
Result IntegrationTestPoints_Set(const TIntegrationFn& IntegrationFn, float actualValue, const char* name, const TPointFn& PointFn, size_t pointCount_)
{
	printf("  %s\n", name);

	// prepare our results
	Result ret;
	ret.name = name;
	ret.values.resize(pointCount_);

	for (size_t pointCount = 1; pointCount <= pointCount_; pointCount++)
	{
		// Generate points
		std::vector<float> points;
		PointFn(points, pointCount);

		// Integerate
		float avgY = 0.0f;
		for (size_t index = 0; index < pointCount; ++index)
		{
			float y = IntegrationFn(points[index]);
			avgY = Lerp(avgY, y, 1.0f / float(index + 1));
		}

		float error = avgY - actualValue;
		ret.values[pointCount-1] = std::abs(error);
	}

	// return result
	return ret;
}

template <typename TIntegrationFn>
void IntegrationTest(const char* name, const TIntegrationFn& IntegrationFn, float actualValue, size_t pointCount)
{
	printf("%s...\n", name);

	// gather results
	std::vector<Result> results;
	results.push_back(IntegrationTestPoints_Sequence(IntegrationFn, actualValue, "Kritzinger", KritzingerLDS, pointCount));
	results.push_back(IntegrationTestPoints_Sequence(IntegrationFn, actualValue, "GoldenRatio", GoldenRatioLDS, pointCount));
	results.push_back(IntegrationTestPoints_Set(IntegrationFn, actualValue, "Regular", RegularPoints, pointCount));
	results.push_back(IntegrationTestPoints_Set(IntegrationFn, actualValue, "RegularWalls", RegularPointsTouchWalls, pointCount));
	results.push_back(IntegrationTestPoints_Set(IntegrationFn, actualValue, "Stratified", Stratified, pointCount));

	// Write CSV
	char fileName[1024];
	sprintf_s(fileName, "out/%s.csv", name);
	FILE* file = nullptr;
	fopen_s(&file, fileName, "wb");
	if (file)
	{
		// write header
		fprintf(file, "\"Samples\"");
		for (const Result& result : results)
			fprintf(file, ",\"%s\"", result.name);
		fprintf(file, "\n");

		// write the data
		for (size_t i = 0; i < results[0].values.size(); ++i)
		{
			fprintf(file, "\"%i\"", (int)i);
			for (const Result& result : results)
				fprintf(file, ",\"%f\"", result.values[i]);
			fprintf(file, "\n");
		}

		fclose(file);
	}
}

int main(int argc, char** argv)
{
	_mkdir("out");

	static const int c_pointCount = 256;

	IntegrationTest("Triangle", Function_Triangle, 0.5f, c_pointCount);
	IntegrationTest("Step", Function_Step, 0.6f, c_pointCount);
	IntegrationTest("Sine", Function_Sine, 2.0f / c_pi, c_pointCount);

	return 0;
}
/*
TODO:
- could try and do what the paper did for the quadratic terms thing, so it is N^2, not N^3 for generating points.
 * or just mention it
- make numberline images
 * nah?
*/