from fastapi import FastAPI, File, UploadFile, BackgroundTasks
from uuid import uuid4
from app.storage import upload_files_to_s3
from app.database import save_pipeline_run, get_pipeline_run_status
from fastapi.middleware.cors import CORSMiddleware
import subprocess  # Import subprocess for running Nextflow
import os
import boto3
import json

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.post("/upload")
async def upload_fastq(
    background_tasks: BackgroundTasks,
    files: list[UploadFile] = File(...)
):
    run_id = str(uuid4())
    s3_paths = upload_files_to_s3(run_id, files)
    save_pipeline_run(run_id, status="running", results=None)
    background_tasks.add_task(trigger_pipeline, run_id, s3_paths)
    return {"run_id": run_id, "status": "Pipeline started"}

def trigger_pipeline(run_id: str, s3_paths: list[str]):
    try:
        print(f"Starting EC2 launch for run_id: {run_id}")

        ec2 = boto3.client('ec2')

        params_json = json.dumps({"s3_paths": s3_paths})

        user_data = f"""#!/bin/bash
        mkdir /home/ubuntu/runs/{run_id}
        echo '{params_json}' > /home/ubuntu/runs/{run_id}/params.json
        nextflow run /path/to/main.nf -params-file /home/ubuntu/runs/{run_id}/params.json -work-dir /home/ubuntu/runs/{run_id}/work
        aws s3 sync /home/ubuntu/runs/{run_id}/work/results/ s3://your-results-bucket/{run_id}/results/
        """

        response = ec2.run_instances(
            LaunchTemplate={'LaunchTemplateId': 'lt-02f0dfa62a554010a'},
            MinCount=1,
            MaxCount=1,
            UserData=user_data
        )

        print(f"AWS EC2 API Response: {response}") # Added print
        instance_id = response['Instances'][0]['InstanceId']
        print(f"EC2 instance launched: {instance_id}")

    except Exception as e:
        print(f"Error launching EC2 instance: {e}")
        # save_pipeline_run(run_id, status="failed", results=str(e)) # commented out

#@app.get("/results/{run_id}")
#async def get_results(run_id: str):
#    status = get_pipeline_run_status(run_id)
#    return status
#
#@app.post("/update_pipeline_status")
#async def update_pipeline_status(run_data: dict):
#    run_id = run_data.get("run_id")
#    status = run_data.get("status")
#    if run_id and status:
#        save_pipeline_run(run_id, status=status, results="...") # add results here if needed
#        return {"message": "Pipeline status updated"}
#    else:
#        return {"message": "Invalid request"}, 400