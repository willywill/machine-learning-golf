using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[ExecuteInEditMode]
public class RotateBlades : MonoBehaviour
{
    public float speed = 0.0f;

    // Update is called once per frame
    private void Update()
    {
        transform.eulerAngles = new Vector3(transform.eulerAngles.x, transform.eulerAngles.y, transform.eulerAngles.z + (Time.deltaTime * speed));
    }
}
