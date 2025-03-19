from fastapi import FastAPI, File, UploadFile, BackgroundTasks
from uuid import uuid4
from app.storage import upload_files_to_s3
from app.pipeline import trigger_pipeline
from app.database import save_pipeline_run, get_pipeline_run_status
from fastapi.middleware.cors import CORSMiddleware  # Import CORS middleware

app = FastAPI()

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],  # Replace with your frontend's origin
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
    
    # Upload files to S3 ## Modify
    s3_paths = upload_files_to_s3(run_id, files)
    
    # Save initial state to DB
    save_pipeline_run(run_id, status="running", results=None)
    
    # Trigger pipeline asynchronously ## Modify
    #background_tasks.add_task(trigger_pipeline, run_id, s3_paths)
    
    return {"run_id": run_id, "status": "Pipeline started"}

@app.get("/results/{run_id}")
async def get_results(run_id: str):
    status = get_pipeline_run_status(run_id)
    return status