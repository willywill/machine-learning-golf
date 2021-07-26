using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// TODO: Make line length tied to the balls putting power
// TODO: Make a vertical power bar GUI
// TODO: Introduce environmental factos - wind, gravity?
// TODO: Add a GUI for the number of putts
// TODO: Move these todo's to the readme, or github issues, or something...
public class BallController : MonoBehaviour
{
    public float changeAngleSpeed = 40.0f;
    public float lineLength = 0.5f;
    public float maxPower = 0.5f;

    public Rigidbody ball;
    public LineRenderer line;
    private float angle = 180.0f;

    private void Start()
    {
        ball = GetComponent<Rigidbody>();
        ball.maxAngularVelocity = 1000;

        line = GetComponent<LineRenderer>();
    }

    // Update is called once per frame
    private void Update()
    {
        if (Input.GetKey(KeyCode.LeftArrow)) 
        {
            ChangeAngle(-1);
        }

        if (Input.GetKey(KeyCode.RightArrow))
        {
            ChangeAngle(1);
        }

        if (Input.GetKey(KeyCode.Space))
        {
            PuttBall();
        }

        SetProjectedPath();
    }

    private void ChangeAngle(int direction)
    {
        angle += changeAngleSpeed * Time.deltaTime * direction;
    }

    private void SetProjectedPath() 
    {
        line.SetPosition(0, transform.position);
        line.SetPosition(1, transform.position + Quaternion.Euler(0, angle, 0) * Vector3.forward * lineLength);
    }

    private void PuttBall()
    {
        ball.AddForce(Quaternion.Euler(0, angle, 0) * Vector3.forward * maxPower, ForceMode.Impulse);
    }
}
