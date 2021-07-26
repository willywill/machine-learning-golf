using System.Threading.Tasks;
using UnityEngine;

// [ExecuteInEditMode]
public class LowPolyClouds : MonoBehaviour
{
    private struct Cloud
    {
        public float Scale;
        public int X;
        public int Y;
    }

    private struct Batch
    {
        public int Length;
        public Cloud[] Clouds;
        public Matrix4x4[] Objects;
    }

    private struct FrameParams
    {
        public float Time;
        public float DeltaTime;
        public int BatchIndex;
    }

    private void Start()
    {
        // total count of clouds
        int count = m_cloudsCount * m_cloudsCount;

        // total count of batches. each size limited by unity to 1023
        int batchCount = count / BatchSize + 1;

        // initialize data arrays
        m_batches = new Batch[batchCount];
        m_tasks = new Task[batchCount];

        for (int i = 0; i < batchCount; i++)
        {
            // last array lenght can be less than 1023
            int length = Mathf.Min(BatchSize, count - i * BatchSize);
            m_batches[i].Length = length;
            m_batches[i].Clouds = new Cloud[length];
            m_batches[i].Objects = new Matrix4x4[length];
        }

        // pivot of clouds should be at center. so just shift each cloud 
        var offset = -m_cloudsCount * 0.5f;
        // initialize data for each cloud
        for (int cloudY = 0; cloudY < m_cloudsCount; cloudY++)
        {
            for (int cloudX = 0; cloudX < m_cloudsCount; cloudX++)
            {
                var cloud = new Cloud
                {
                    Scale = 0,
                    X = cloudX,
                    Y = cloudY,
                };

                // position of this cloud in world
                var position = new Vector3
                {
                    x = offset + transform.position.x + cloudX * m_cloudSize,
                    y = transform.position.y,
                    z = offset + transform.position.z + cloudY * m_cloudSize,
                };

                // convert X and Y to batch indices
                int index = cloudY * m_cloudsCount + cloudX;

                int x = index / BatchSize;
                int y = index % BatchSize;

                m_batches[x].Clouds[y] = cloud;
                m_batches[x].Objects[y] = Matrix4x4.TRS(position, Quaternion.identity, Vector3.zero);
            }
        }
    }

    private void Update()
    {
        // each batch will be updated in separate thread
        for (int batchIndex = 0; batchIndex < m_batches.Length; batchIndex++)
        {
            // to avoid allocations while creating delegates like ()=> UpdateBatch(frameParams)
            // i'm sending frame params as object and then cast it to FrameParams
            FrameParams frameParams = new FrameParams
            {
                BatchIndex = batchIndex,
                DeltaTime = Time.deltaTime,
                Time = Time.time
            };

            // creation of new task allocates some memory but i can't do anythung with it
            m_tasks[batchIndex] = Task.Factory.StartNew(UpdateBatch, frameParams);
        }

        // just wait all tasks to be completed
        Task.WaitAll(m_tasks);

        // and send all objects to render
        for (int batchIndex = 0; batchIndex < m_batches.Length; batchIndex++)
        {
            Graphics.DrawMeshInstanced(m_mesh, 0, m_material, m_batches[batchIndex].Objects);
        }
    }

    private void UpdateBatch(object input)
    {
        FrameParams frameParams = (FrameParams)input;
        for (int cloudIndex = 0; cloudIndex < m_batches[frameParams.BatchIndex].Length; cloudIndex++)
        {
            // just to simplify code below
            int i = frameParams.BatchIndex; int j = cloudIndex;

            // calculate noise based on coordinates of cloud and current time
            float x = m_batches[i].Clouds[j].X * m_texScale + frameParams.Time * m_timeScale;
            float y = m_batches[i].Clouds[j].Y * m_texScale + frameParams.Time * m_timeScale;
            float noise = Mathf.PerlinNoise(x, y);

            // based on noise we can understand scale direction
            int dir = noise > m_minNoiseSize ? 1 : -1;

            // calcuate new scale and clamp it
            float shift = m_scaleSize * frameParams.DeltaTime * dir;
            float scale = m_batches[i].Clouds[j].Scale + shift;
            scale = Mathf.Clamp(scale, 0, m_maxScale);
            m_batches[i].Clouds[j].Scale = scale;

            // set new scale to object matrix
            m_batches[i].Objects[j].m00 = scale;
            m_batches[i].Objects[j].m11 = scale;
            m_batches[i].Objects[j].m22 = scale;
        }
    }

    [SerializeField]
    private Mesh m_mesh;
    [SerializeField]
    private Material m_material;
    [SerializeField]
    private float m_cloudSize = 1;
    [SerializeField]
    private float m_maxScale = 1;
    [SerializeField]
    private float m_timeScale = 0.05f;
    [SerializeField]
    private float m_texScale = 0.1f;
    [SerializeField]
    private float m_minNoiseSize = 0.6f;
    [SerializeField]
    private float m_scaleSize = 1.5f;
    [SerializeField]
    private int m_cloudsCount = 100;

    // maximum size of matrices array that can be passed to Graphics.DrawMeshInstanced
    private const int BatchSize = 1023;

    private Batch[] m_batches;
    private Task[] m_tasks;
}