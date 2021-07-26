Shader "iCE/LowPolyClouds"
{
    Properties
    {
        // _Color ("Color", Color) = (1,1,1,1)
        _MieColor ("Mie Color", Color) = (1,1,1,1)
        _SkyColor ("Sky Color", Color) = (1,1,1,1)
        _GroundColor ("Ground Color", Color) = (1,1,1,1)
        _G ("Mie Anisotrophy", Range(0,1)) = 0.5
        _MieScatt ("Mie Scattering", Range(0,10)) = 2.5
        _AmbientInt ("Ambient Intensity", Range(0,10)) = 1.0
        _Extinct ("Extinction", Range(0,10)) = 1.0
    }
    SubShader
    {
        Tags {"Queue" = "Transparent" "RenderType"="Transparent" "ForceNoShadowCasting" = "True" }
        LOD 200

        CGPROGRAM
        // Physically based Standard lighting model, and enable shadows on all light types
        #pragma surface surf Standard fullforwardshadows alpha:fade

        // Use shader model 3.0 target, to get nicer looking lighting
        #pragma target 3.0

        sampler2D _MainTex;

        struct Input
        {
            float3 viewDir;
            float3 worldNormal;
            float4 screenPos;
            float eyeDepth;
        };

        half _G;
        half _MieScatt;
        half _AmbientInt;
        half _Extinct;
        fixed4 _MieColor;
        fixed4 _GroundColor;
        fixed4 _SkyColor;
        sampler2D_float _CameraDepthTexture;
        float4 _CameraDepthTexture_TexelSize;

        /*
        MiePhase
        */
        float MiePhase(float cosTheta, float g)
        {
            return (1.0 / (4.0 * UNITY_PI)) * ((1.0 - pow(g, 2)) / pow(1.0 - 2.0 * g * cosTheta + pow(g, 2), 1.5));
        }

        float clamp01(float value) 
        {
            return clamp(value, 0.0, 1.0);
        }

        // Add instancing support for this shader. You need to check 'Enable Instancing' on materials that use the shader.
        // See https://docs.unity3d.com/Manual/GPUInstancing.html for more information about instancing.
        // #pragma instancing_options assumeuniformscaling
        UNITY_INSTANCING_BUFFER_START(Props)
            // put more per-instance properties here
        UNITY_INSTANCING_BUFFER_END(Props)

        void vert (inout appdata_full v, out Input o)
        {
            UNITY_INITIALIZE_OUTPUT(Input, o);
            COMPUTE_EYEDEPTH(o.eyeDepth);
        }

        void surf (Input IN, inout SurfaceOutputStandard o)
        {
            float3 up = float3(0.0, 1.0, 0.0);
            float3 rayDir = normalize(IN.viewDir);
            float3 lightDir = normalize(_WorldSpaceLightPos0.xyz);

            float3 skyCol = (dot(IN.worldNormal, up) * 0.5 + 0.5) * _SkyColor.rgb;
            float3 groundCol = (dot(IN.worldNormal, -up) * 0.5 + 0.5) * _GroundColor.rgb;
            float3 mieCol = (MiePhase(dot(rayDir, -lightDir), _G) + MiePhase(dot(rayDir, lightDir), max(0.1, (1.0 - _G)))) * _MieScatt;
            mieCol *= _MieColor.rgb;

            fixed4 ambientCol = fixed4(skyCol + groundCol, 0.0) * _AmbientInt;
            o.Albedo = ambientCol + mieCol;

            o.Metallic = 0.0;
            o.Smoothness = 0.0;

            float rawZ = SAMPLE_DEPTH_TEXTURE_PROJ(_CameraDepthTexture, UNITY_PROJ_COORD(IN.screenPos));
            float depth = Linear01Depth(rawZ);

            o.Alpha = exp(-_Extinct * depth);
        }
        ENDCG
    }
    FallBack "Diffuse"
}
