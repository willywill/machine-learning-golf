#ifndef ATMOSPHERIC_SCATTERING_INCLUDED
#define ATMOSPHERIC_SCATTERING_INCLUDED

/* Definitions */
#define mieExtinctionFactor 1.11
#define mieScatteringCoeff 21e-6
#define rayleighExtinction float3(5.8e-6, 13.5e-6, 33.1e-6)
#define ozoneExtinction float3(3.426e-7, 8.298e-7, 0.356e-7)
#define starHash float4(641, -113, 271, 1117)

/* Constants */
// uniform float4 ScatteringProfile_0;
// uniform float4 ScatteringProfile_1;
// uniform float4 ScatteringProfile_2;

// #define _PlanetRadius ScatteringProfile_0.x
// #define _AtmosphereRadius ScatteringProfile_0.y
// #define _RayleighScattering ScatteringProfile_0.z
// #define _MieScattering ScatteringProfile_0.w

// #define _G ScatteringProfile_1.x
// #define _OzoneExtinctionFactor ScatteringProfile_1.y
// #define _SunIlluminance ScatteringProfile_1.z
// #define _SkyIlluminance ScatteringProfile_1.w

// #define _SunSize ScatteringProfile_2.x
// #define _StarSize ScatteringProfile_2.y
// #define _StarDensity ScatteringProfile_2.z
// #define _LightPollutionFactor ScatteringProfile_2.w

uniform float _PlanetRadius;
uniform float _AtmosphereRadius;
uniform float _RayleighScattering;
uniform float _MieScattering;
uniform float _G;
uniform float _OzoneExtinctionFactor;
uniform float _SunIlluminance;
uniform float _SkyIlluminance;
uniform float _SunSize;
uniform float _StarSize;
uniform float _StarDensity;
uniform float _LightPollutionFactor;

uniform float _DistanceScale;
uniform float _HeightFog;

uniform float4 _CameraPos;
// uniform float4 _LightDir;

// Frustrum Corners
uniform float4 _TopLeftCorner;
uniform float4 _TopRightCorner;
uniform float4 _BottomLeftCorner;
uniform float4 _BottomRightCorner;

uniform sampler2D _TransmittanceLUT;
uniform sampler3D _SkyboxRayleighLUT;
uniform sampler3D _SkyboxMieLUT;
uniform sampler3D _Inscattering;
uniform sampler3D _Extinction;

RWTexture3D<float4> _InscatteringLUT;
RWTexture3D<float4> _ExtinctionLUT;

/* Samples */
#define primaryRaySteps 16
#define secondaryRaySteps 8
#define primaryRayStepGrowth 1.15
#define secondaryRayStepGrowth 1.15

/* Structs */
struct SkyParams
{
    float planetRadius;
    float atmosphereRadius;
    float rayleighScattering;
    float mieScattering;
    float mieAnisotrophy;
    float ozoneLayerDensity;
    float sunIlluminance;
    float skyIlluminance;
    float sunSize;
    float starSize;
    float starDensity;
    float lightPollutionFactor;
};

struct SkyContext
{
    float cosTheta;
    float cosBeta;
    float cosAlpha;
    float3 rayDir;
    float3 rayOrigin;
    float3 sunDir;
};

struct SkyLutContext {
    float3 rayDir;
    float3 lightDir;
    float height;
};

struct Scattering {
    float3 inscatter;
    float3 extinction;
};

SkyLutContext EncodeSkyboxLUTCoords(SkyParams skyParams, uint3 id, float w, float h, float d)
{
    SkyLutContext skyLutContext;

    float3 coords = float3(id.x / (w - 1), id.y / (h - 1), id.z / (d - 1));

    float atmosphereHeight = skyParams.atmosphereRadius - skyParams.planetRadius;

    float height = coords.x * coords.x * atmosphereHeight;
    float cosHorizon = -(sqrt(height * (2 * skyParams.planetRadius + height)) / (skyParams.planetRadius + height));

    float viewZenithAngle = coords.y;

    if (viewZenithAngle > 0.5)
    {
        viewZenithAngle = cosHorizon + pow(abs(viewZenithAngle - 0.5) * 2, 5) * (1 - cosHorizon);
    }
    else
    {
        viewZenithAngle = cosHorizon - pow(abs(viewZenithAngle * 2), 5) * (1 + cosHorizon);
    }

    float sunZenithAngle = (tan((2 * coords.z - 1 + 0.26) * 0.75)) / (tan(1.26 * 0.75));

    float3 rayDir = float3(sqrt(saturate(1 - viewZenithAngle * viewZenithAngle)), viewZenithAngle, 0);
    float3 lightDir = -float3(sqrt(saturate(1 - sunZenithAngle * sunZenithAngle)), sunZenithAngle, 0);

    skyLutContext.rayDir = rayDir;
    skyLutContext.lightDir = lightDir;
    skyLutContext.height = height + skyParams.planetRadius; // camera height between ground level 0 and atsmosphere height

    return skyLutContext;
}

float3 DecodeSkyboxLUTCoords(SkyParams skyParams, float3 rayDir, float3 lightDir, float3 cameraPos)
{
    float3 rayOrigin = float3(0.0, max(0.001, cameraPos.y), 0.0);
    float atmosphereHeight = skyParams.atmosphereRadius - skyParams.planetRadius;

    float height = length(rayOrigin);
    float3 normal = normalize(rayOrigin);

    float viewZenith = dot(normal, rayDir);
    float sunZenith = dot(normal, -lightDir);

    float3 coords = float3(height / atmosphereHeight, viewZenith * 0.5 + 0.5, sunZenith * 0.5 + 0.5);

    coords.x = pow(height / atmosphereHeight, 0.5);
    float cosHorizon = -(sqrt(height * (2 * skyParams.planetRadius + height)) / (skyParams.planetRadius + height));

    if (viewZenith > cosHorizon)
    {
        coords.y = 0.5 * pow((viewZenith - cosHorizon) / (1 - cosHorizon), 0.2) + 0.5;
    }
    else
    {
        coords.y = 0.5 * pow((cosHorizon - viewZenith) / (cosHorizon + 1), 0.2);
    }

    coords.z = 0.5 * ((atan(max(sunZenith, -0.1975) * tan(1.26 * 1.1)) / 1.1) + (1 - 0.26));

    return coords;
}

/*
  RayleighPhase
*/
float RayleighPhase(float cosTheta) 
{
    return (3.0 / (16.0 * UNITY_PI)) * (1.0 + (cosTheta * cosTheta));
}

/*
  MiePhase
*/
float MiePhase(float cosTheta, float g)
{
    return (1.0 / (4.0 * UNITY_PI)) * ((1.0 - pow(g, 2)) / pow(1.0 - 2.0 * g * cosTheta + pow(g, 2), 1.5));
}

/*
  SunDisk
  TODO: Use physically based approach - https://media.contentapi.ea.com/content/dam/eacom/frostbite/files/s2016-pbs-frostbite-sky-clouds-new.pdf
*/
float3 SunDisk(SkyContext skyContext, SkyParams skyParams, float3 transmittance, float limbDarkeningFactor)
{
    float3 delta = skyContext.sunDir - skyContext.rayDir;
    float dist = length(delta);
    float softenRation = skyParams.sunSize * limbDarkeningFactor;
    float limbDarkening = 1.0 - smoothstep(softenRation, skyParams.sunSize, dist);

    return (skyParams.sunIlluminance * limbDarkening) * transmittance;
}

/*
  StarPattern - based on common random hash
*/
float StarPattern(int n)
{
    return frac(sin(float(n) * 543.21) * 43758.5453);
}

/*
  IntegrateBlackbody
*/
float IntegrateBlackbody(float x) 
{ 
    return (6.0 + x * (6.0 + x * (3.0 + x))) * exp(-x);
}

/*
  PlankBlackbody
*/
float PlankBlackbody(float T, float lambda1, float lambda0) 
{
	  float A = 1.1;
    float B = 1.0 / 1.05;
	  float C0 = 0.014387770;
    float C = C0 / (B * T);

	return 100.0 * A / B * pow(100.0 / C0, 4.0) * (IntegrateBlackbody(C / lambda1) - IntegrateBlackbody(C / lambda0));
}

/*
  StarColor - https://www.shadertoy.com/view/XdsGWs
*/
float3 StarColor(float temp)
{
    float r = PlankBlackbody(temp, 7e-6, .55e-6);
    float g = PlankBlackbody(temp, .55e-6, .49e-6);
    float b = PlankBlackbody(temp, .49e-6, .4e-6);

    return float3(r, g, b) * 1e-14;
}

/*
  ProceduralStars
*/
float3 ProceduralStars(SkyContext skyContext, SkyParams skyParams)
{
    float3 starPos = skyContext.rayDir / skyParams.starSize;
    // Move stars along as the moon moves.
    // starPos.z += skyContext.sunDir.y * 10;
    // starPos.y += skyContext.sunDir.y * 10;
    float3 center = round(starPos);
    float hash = dot(starHash.xyz, center) % starHash.w;
    float threshold = starHash.w * skyParams.starDensity;

    if (abs(hash) < threshold)
    {
      float dist = length(starPos - center);
      float3 color = StarColor(40000.0 * exp(-3.0 * StarPattern(6 * hash + 5)));

      return max(0.0, saturate(color * pow(saturate(0.5 - dist * dist) * 2, 14)));
    }

    return float3(0.0, 0.0, 0.0);
}

/*
  RaySphereIntersection
*/
float2 RaySphereIntersection(float3 rayOrigin, float3 rayDir, float sphereRadius)
{
    float a = dot(rayDir, rayDir);
    float b = 2.0 * dot(rayOrigin, rayDir);
    float c = dot(rayOrigin, rayOrigin) - (sphereRadius * sphereRadius);
    float d = b * b - 4 * a * c;
    if (d < 0)
    {
        return float2(1e5, -1e5);
    }
    else
    {
        d = sqrt(d);
        return float2(-b - d, -b + d) / (2 * a);
    }
}

/*
  GeometricSeries
  Sum of first n terms in a geometric progression starting with 1
  a * ( 1 - r^n ) / ( 1 - r ), here a is 1.
*/
float GeometricSeries(float commonRatio, float numTerms) 
{
    return (1.0 - pow(commonRatio, numTerms)) / (1.0 - commonRatio);
}

/*
  ComputeTransmittance
*/
float3 ComputeTransmittance(SkyParams skyParams, float3 primaryOpticalDepthRayleigh, float3 primaryOpticalDepthMie, float3 secondaryOpticalDepthRayleigh, float3 secondaryOpticalDepthMie) 
{
    float3 density = float3(0.0, 0.0, 0.0);

    float3 mieDensity = mieScatteringCoeff * (primaryOpticalDepthMie + secondaryOpticalDepthMie) * mieExtinctionFactor;
    float3 rayleighDensity = rayleighExtinction * (primaryOpticalDepthRayleigh + secondaryOpticalDepthRayleigh);
    float3 ozoneDensity = ozoneExtinction * (primaryOpticalDepthRayleigh + secondaryOpticalDepthRayleigh) * skyParams.ozoneLayerDensity;

    density = mieDensity + (rayleighDensity + ozoneDensity);

    return exp(-density);
}

/*
  ComputeStars
*/
float3 ComputeStars(SkyParams skyParams, SkyContext skyContext) 
{
    // float lightPollution = lerp(0.0, 1.0, rayleighScaleHeight / 350000);
    float starPossibility = lerp(0.0, skyContext.cosAlpha, max(0.0, -skyContext.sunDir.y));
    float3 stars = ProceduralStars(skyContext, skyParams) * starPossibility;
    return stars;
}

void ApplyHeightFog(float3 worldPos, inout float density)
{
    density *= saturate(exp(-(worldPos.y + 0.0) * _HeightFog));
}

/*
  ComputeAerialPerspective
*/
float4 ComputeAerialPerspective(SkyParams skyParams, float3 coords, float4 color)
{
    float4 inscattering = tex3D(_Inscattering, coords);
    float4 extinction = tex3D(_Extinction, coords);

    if (coords.z > 0.99999)
    {
        inscattering = 0;
        extinction = 1;
        // color = 0; // For debugging
    }

    return color * extinction + (inscattering * skyParams.skyIlluminance * 2);
}

/*
  ApplyScaleHeight
*/
void ApplyScaleHeight(SkyContext skyContext, SkyParams skyParams, out float rayleighScaleHeight, out float mieScaleHeight) 
{
    rayleighScaleHeight = skyParams.rayleighScattering * lerp(1.0, skyParams.lightPollutionFactor, max(0.0, skyContext.sunDir.y));
    mieScaleHeight = skyParams.mieScattering * lerp(1.0, skyParams.lightPollutionFactor, max(0.0, skyContext.sunDir.y));
}

/*
  IntegrateScattering
*/
void IntegrateScattering(SkyContext skyContext, SkyParams skyParams, out float3 totalRayleigh, out float3 totalMie, out float3 transmittance)
{
    // Calculate the step size of the primary ray.
    float2 position = RaySphereIntersection(skyContext.rayOrigin, skyContext.rayDir, skyParams.atmosphereRadius);

    position.y = min(position.y, RaySphereIntersection(skyContext.rayOrigin, skyContext.rayDir, skyParams.planetRadius).x);

    float dist = position.y - position.x;
    float stepSize = dist / GeometricSeries(primaryRayStepGrowth, float(primaryRaySteps));

    float mieScaleHeight = 0.0;
    float rayleighScaleHeight = 0.0;
    float primaryOpticalDepthRayleigh = 0.0;
    float primaryOpticalDepthMie = 0.0;

    ApplyScaleHeight(skyContext, skyParams, rayleighScaleHeight, mieScaleHeight);

    float primaryTime = 0.0;

    // Sample the primary ray.
    [loop]
    for (int i = 0; i < primaryRaySteps; i++) 
    {
      // Calculate the primary ray sample position.
      float3 primaryRayPos = skyContext.rayOrigin + skyContext.rayDir * (primaryTime + stepSize * 0.5);

      // Calculate the height of the sample.
      float primaryHeight = length(primaryRayPos) - skyParams.planetRadius;

      // Calculate the optical depth of the Rayleigh and Mie scattering for this step.
      float opticalDepthRayleighStep = exp(-primaryHeight / rayleighScaleHeight) * stepSize;
      float opticalDepthMieStep = exp(-primaryHeight / mieScaleHeight) * stepSize;

      // Accumulate optical depth.
      primaryOpticalDepthRayleigh += opticalDepthRayleighStep;
      primaryOpticalDepthMie += opticalDepthMieStep;

      // Calculate the step size of the secondary ray.
      float secondaryStepSize = RaySphereIntersection(primaryRayPos, skyContext.sunDir, skyParams.atmosphereRadius).y / GeometricSeries(secondaryRayStepGrowth, float(secondaryRaySteps));

      // Initialize the secondary ray time.
      float secondaryTime = 0.0;

      // Initialize optical depth accumulators for the secondary ray.
      float secondaryOpticalDepthRayleigh = 0.0;
      float secondaryOpticalDepthMie = 0.0;

      // Sample the secondary ray.
      [loop]
      for (int j = 0; j < secondaryRaySteps; j++) 
      {
        // Calculate the secondary ray sample position.
        float3 secondaryRayPos = primaryRayPos + skyContext.sunDir * (secondaryTime + secondaryStepSize * 0.5);

        // Calculate the height of the sample.
        float secondaryHeight = length(secondaryRayPos) - skyParams.planetRadius;

        // Accumulate the optical depth.
        secondaryOpticalDepthRayleigh += exp(-secondaryHeight / rayleighScaleHeight) * secondaryStepSize;
        secondaryOpticalDepthMie += exp(-secondaryHeight / mieScaleHeight) * secondaryStepSize;

        // Increment the secondary ray time.
        secondaryTime += secondaryStepSize;
        secondaryStepSize *= secondaryRayStepGrowth;
      }

      // Calculate attenuation.
      transmittance = ComputeTransmittance(skyParams, primaryOpticalDepthRayleigh, primaryOpticalDepthMie, secondaryOpticalDepthRayleigh, secondaryOpticalDepthMie);

      // Accumulate scattering.
      totalRayleigh += opticalDepthRayleighStep * transmittance;
      totalMie += opticalDepthMieStep * transmittance;

      // Increment the primary ray time.
      primaryTime += stepSize;
      stepSize *= primaryRayStepGrowth;
    }
}

/*
  IntegrateScatteringPrecomputed
  TODO: Use analytical scattering term from DICE?
*/
void IntegrateScatteringPrecomputed(SkyContext skyContext, SkyParams skyParams, uint3 coords, uint sampleCount)
{
    // Calculate the step size of the primary ray.
    float2 position = RaySphereIntersection(skyContext.rayOrigin, skyContext.rayDir, skyParams.atmosphereRadius);

    position.y = min(position.y, RaySphereIntersection(skyContext.rayOrigin, skyContext.rayDir, skyParams.planetRadius).x);

    float dist = (position.y - position.x) * 0.5;
    float stepSize = dist / GeometricSeries(primaryRayStepGrowth, float(sampleCount));
    stepSize *= _DistanceScale;

    float mieScaleHeight = 0.0;
    float rayleighScaleHeight = 0.0;
    float primaryOpticalDepthRayleigh = 0.0;
    float primaryOpticalDepthMie = 0.0;

    float3 totalRayleigh = float3(0.0, 0.0, 0.0);
    float3 totalMie = float3(0.0, 0.0, 0.0);
    float3 transmittance = float3(0.0, 0.0, 0.0);

    _InscatteringLUT[coords] = float4(0, 0, 0, 1);
	  _ExtinctionLUT[coords] = float4(1, 1, 1, 1);

    ApplyScaleHeight(skyContext, skyParams, rayleighScaleHeight, mieScaleHeight);

    float primaryTime = 0.0;

    // Sample the primary ray
    [loop]
    for (coords.z = 1; coords.z < sampleCount; coords.z += 1)
    {
      // Calculate the primary ray sample position.
      float3 primaryRayPos = skyContext.rayOrigin + skyContext.rayDir * (primaryTime + stepSize * 0.5);

      // Calculate the height of the sample.
      float primaryHeight = length(primaryRayPos) - skyParams.planetRadius;

      // Calculate the optical depth of the Rayleigh and Mie scattering for this step.
      float opticalDepthRayleighStep = exp(-primaryHeight / rayleighScaleHeight) * stepSize;
      float opticalDepthMieStep = exp(-primaryHeight / mieScaleHeight) * stepSize;

      // Accumulate optical depth.
      primaryOpticalDepthRayleigh += opticalDepthRayleighStep;
      primaryOpticalDepthMie += opticalDepthMieStep;

      // Calculate the step size of the secondary ray.
      float secondaryStepSize = RaySphereIntersection(primaryRayPos, skyContext.sunDir, skyParams.atmosphereRadius).y / GeometricSeries(secondaryRayStepGrowth, float(sampleCount / 2));

      // Initialize the secondary ray time.
      float secondaryTime = 0.0;

      // Initialize optical depth accumulators for the secondary ray.
      float secondaryOpticalDepthRayleigh = 0.0;
      float secondaryOpticalDepthMie = 0.0;

      // Sample the secondary ray.
      [loop]
      for (int j = 0; j < (sampleCount / 2); j++) 
      {
        // Calculate the secondary ray sample position.
        float3 secondaryRayPos = primaryRayPos + skyContext.sunDir * (secondaryTime + secondaryStepSize * 0.5);

        // Calculate the height of the sample.
        float secondaryHeight = length(secondaryRayPos) - skyParams.planetRadius;

        // Accumulate the optical depth.
        secondaryOpticalDepthRayleigh += exp(-secondaryHeight / rayleighScaleHeight) * secondaryStepSize;
        secondaryOpticalDepthMie += exp(-secondaryHeight / mieScaleHeight) * secondaryStepSize;

        // Increment the secondary ray time.
        secondaryTime += secondaryStepSize;
        secondaryStepSize *= secondaryRayStepGrowth;
      }

      // Calculate attenuation.
      transmittance = ComputeTransmittance(skyParams, primaryOpticalDepthRayleigh, primaryOpticalDepthMie, secondaryOpticalDepthRayleigh * 0.1, secondaryOpticalDepthMie);

      // Accumulate scattering.
      totalRayleigh += opticalDepthRayleighStep * transmittance;
      totalMie += opticalDepthMieStep * transmittance;

      // Increment the primary ray time.
      primaryTime += stepSize;
      stepSize *= primaryRayStepGrowth;

      float phaseRayleigh = max(0.0, RayleighPhase(skyContext.cosTheta));
      float phaseMie = max(0.0, MiePhase(skyContext.cosTheta, skyParams.mieAnisotrophy));

      float3 rayleighTerm = phaseRayleigh * rayleighExtinction * totalRayleigh * skyParams.skyIlluminance;
      float3 mieTerm = phaseMie * mieScatteringCoeff * totalMie;
      float3 inscattering = rayleighTerm + mieTerm;

      inscattering = max(float3(0.0, 0.0, 0.0), inscattering);
      transmittance =  max(float3(0.0, 0.0, 0.0), transmittance);

      // Compute inscattering and extinction for each slice in the 3D texture
      _InscatteringLUT[coords] = float4(inscattering, 1);
		  _ExtinctionLUT[coords] = float4(transmittance, 1);
    }
}

/*
  PrecomputeDirectLight
*/
float4 PrecomputeDirectLight(SkyContext skyContext, SkyParams skyParams)
{
    float3 totalRayleigh = float3(0.0, 0.0, 0.0);
    float3 totalMie = float3(0.0, 0.0, 0.0);
    float3 transmittance = float3(0.0, 0.0, 0.0);

    IntegrateScattering(skyContext, skyParams, totalRayleigh, totalMie, transmittance);

    return float4(transmittance, 1.0);
}

/*
  AtmospericScattering
*/
Scattering AtmospericScattering(SkyContext skyContext, SkyParams skyParams)
{
    Scattering scatter;

    float3 totalRayleigh = float3(0.0, 0.0, 0.0);
    float3 totalMie = float3(0.0, 0.0, 0.0);
    float3 transmittance = float3(0.0, 0.0, 0.0);

    IntegrateScattering(skyContext, skyParams, totalRayleigh, totalMie, transmittance);

    float phaseRayleigh = RayleighPhase(skyContext.cosTheta);
    float phaseMie = MiePhase(skyContext.cosTheta, skyParams.mieAnisotrophy);

    // Calculate and return the final color.
    float3 rayleighTerm = phaseRayleigh * rayleighExtinction * totalRayleigh * skyParams.skyIlluminance;
    float3 mieTerm = phaseMie * mieScatteringCoeff * totalMie;
    float3 inscattering = rayleighTerm + mieTerm;

    scatter.inscatter = inscattering;
    scatter.extinction = transmittance;

    return scatter;
}

Scattering AtmospericScatteringPrecomputed(SkyContext skyContext, SkyParams skyParams, float3 totalRayleigh, float3 totalMie)
{
    Scattering scatter;

    float phaseRayleigh = RayleighPhase(skyContext.cosTheta);
    float phaseMie = MiePhase(skyContext.cosTheta, skyParams.mieAnisotrophy);

    float3 rayleighTerm = phaseRayleigh * rayleighExtinction * totalRayleigh * skyParams.skyIlluminance;
    float3 mieTerm = phaseMie * mieScatteringCoeff * totalMie;
    float3 inscattering = rayleighTerm + mieTerm;

    scatter.inscatter = inscattering;
    scatter.extinction = totalMie;

    return scatter;
}

#endif // ATMOSPHERIC_SCATTERING_INCLUDED