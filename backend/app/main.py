from fastapi import FastAPI, File, UploadFile, BackgroundTasks
from uuid import uuid4
from app.storage import upload_files_to_s3
from app.database import save_pipeline_run, get_pipeline_run_status
from fastapi.middleware.cors import CORSMiddleware
import subprocess  # Import subprocess for running Nextflow
import os

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
        # Create a directory for the run
        run_dir = f"runs/{run_id}"
        os.makedirs(run_dir, exist_ok=True)

        # Create a params.json file for Nextflow
        params_file = os.path.join(run_dir, "params.json")
        with open(params_file, "w") as f:
            import json
            json.dump({"s3_paths": s3_paths}, f)

        # Run Nextflow pipeline
        subprocess.run(
            [
                "nextflow",
                "run",
                "/Users/other/webapp-project/pipeline/main.nf",  # Path to your Nextflow script
                "-params-file",
                params_file,
                "-work-dir",
                os.path.join(run_dir, "work"),
                "-resume" # Consider using -resume if needed.
            ],
            check=True,
            cwd=run_dir, # Run nextflow in the run directory.
            capture_output=True,
            text=True
        )

        # Handle output files (store in S3, update database, etc.)
        # ... (Implement your logic here)

        save_pipeline_run(run_id, status="completed", results="...")  # Update database
    except subprocess.CalledProcessError as e:
        save_pipeline_run(run_id, status="failed", results=str(e.stderr))
    except Exception as e:
        save_pipeline_run(run_id, status="failed", results=str(e))

@app.get("/results/{run_id}")
async def get_results(run_id: str):
    status = get_pipeline_run_status(run_id)
    return status