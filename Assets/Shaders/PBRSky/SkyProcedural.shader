Shader "iCE/Skybox/PhysicallyBasedSky" 
{
    Properties 
    {
        _PlanetRadius ("Planet Radius", Float) = 6368000
        _AtmosphereRadius ("Atmosphere Radius", Float) = 6471000
        _RayleighScattering ("Rayleigh Scattering", Range (1000, 25000)) = 8000.0
        _MieScattering ("Mie Scattering", Range (1200, 15000)) = 2200
        _G ("Mie Scattering Asymmetry", Range (0.01, 0.999)) = 0.85
        _OzoneExtinctionFactor ("Ozone Absorption Factor", Range (0.001, 25.0)) = 3.0
        _SunIlluminance ("Sun Illuminance", Range (0.1, 10.0)) = 1.0
        _SkyIlluminance ("Sky Illuminance", Range (0.1, 10.0)) = 2.0
        _SunSize ("Sun Size", Range (0.001, 0.1)) = 0.03
        _StarSize ("Star Size", Range (0.001, 0.01)) = 0.005
        _StarDensity ("Star Density", Range (0, 0.005)) = 0.002
        _LightPollutionFactor ("Light Pollution Factor", Range (1.0, 15.0)) = 1.0
    }

    SubShader 
    {
        Tags { "Queue"="Background" "RenderType"="Background" "PreviewType" = "Skybox" }
        Cull Off ZWrite Off

        Pass 
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"
            #include "AtmosphericScattering.cginc"

            uniform sampler2D _VolumeClouds;
            uniform float4 _MainTex_TexelSize;

            struct appdata_t
            {
                float4 vertex : POSITION;
                float3 uv : TEXCOORD0;

                UNITY_VERTEX_INPUT_INSTANCE_ID
            };

            struct v2f
            {
                float4 pos : SV_POSITION;
                float3 vertex : TEXCOORD0;
                float3 uv : TEXCOORD1;

                UNITY_VERTEX_OUTPUT_STEREO
            };

            v2f vert(appdata_t v)
            {
                v2f OUT;
                UNITY_SETUP_INSTANCE_ID(v);
                UNITY_INITIALIZE_VERTEX_OUTPUT_STEREO(OUT);
                OUT.pos = UnityObjectToClipPos(v.vertex);
                OUT.vertex = v.vertex;

                // #if UNITY_UV_STARTS_AT_TOP
				// if (_MainTex_TexelSize.y < 0)
				// OUT.uv.y = 1 - OUT.uv.y;
                // #endif

                OUT.uv = v.uv;
                return OUT;
            }

            SkyContext InitializeSkyContext(float3 rayDir, float3 lightDir, float planetRadius)
            {
              SkyContext ctx;

              float3 up = float3(0.0, 1.0, 0.0);

              ctx.cosTheta = dot(rayDir, lightDir);
              ctx.cosBeta = dot(lightDir, up);
              ctx.cosAlpha = dot(rayDir, up);
              ctx.rayDir = rayDir;
              ctx.rayOrigin = float3(0.0, planetRadius, 0.0);
              ctx.sunDir = lightDir;

              return ctx;
            }

            SkyParams InitializeSkyParams()
            {
                SkyParams params;

                params.planetRadius = _PlanetRadius;
                params.atmosphereRadius = _AtmosphereRadius;
                params.rayleighScattering = _RayleighScattering;
                params.mieScattering = _MieScattering;
                params.mieAnisotrophy = _G;
                params.ozoneLayerDensity = _OzoneExtinctionFactor;
                params.sunIlluminance = _SunIlluminance;
                params.skyIlluminance = _SkyIlluminance;
                params.sunSize = _SunSize;
                params.starSize = _StarSize;
                params.starDensity = _StarDensity;
                params.lightPollutionFactor = _LightPollutionFactor;

                return params;
            }

            float3 GetMoonDirection(float3 lightDir)
            {
              float3 moonDir = -lightDir.xyz;

              moonDir.x += 0.4;
              moonDir.y += 0.31;

              return normalize(moonDir);
            }

            float2 SphereToCartesian(float3 dir)
            {
                float2 lonLat = float2(atan2(-dir.x, dir.z), acos(dir.y));
                return lonLat * float2(1 / radians(360.0), 1 / radians(180.0)) + float2(0.5, 0.0);
            }

            float4 frag(v2f i) : SV_Target
            {
                float3 rayDir = normalize(mul((float3x3)unity_ObjectToWorld, i.vertex));
                // float3 coords = normalize(mul((float3x3)unity_ObjectToWorld, i.uv));
                float3 lightDir = normalize(_WorldSpaceLightPos0.xyz);
                float3 moonDir = GetMoonDirection(lightDir);

                SkyParams skyParams = InitializeSkyParams();

                SkyContext skyContext = InitializeSkyContext(rayDir, lightDir, skyParams.planetRadius);
                SkyContext nightSkyContext = InitializeSkyContext(rayDir, moonDir, skyParams.planetRadius);

                Scattering skyScattering = AtmospericScattering(skyContext, skyParams);
                // Scattering nightSkyScattering = AtmospericScattering(nightSkyContext, skyParams);

                float3 atmosphere = (skyScattering.inscatter * 10) + SunDisk(skyContext, skyParams, skyScattering.extinction, 0.79);
                // float3 atmosphereNight = (nightSkyScattering.inscatter * 0.5) + saturate(ComputeStars(nightSkyContext, skyParams)) + saturate(SunDisk(nightSkyContext, skyParams, nightSkyScattering.extinction, 0.79));
                // return float4(SphereToCartesian(i.uv), 0, 0);
                // float2 coords = SphereToCartesian(rayDir);
                // coords.y = 1.0 - coords.y;
                // float4 cloud = saturate(tex2D(_VolumeClouds, coords.xy) * 1);
                // return cloud;
                // atmosphere = atmosphere * (1.0 - cloud.a) + cloud.rgb;

                return float4(atmosphere, 1.0);
            }
            ENDCG
        }
    }

    Fallback Off
}