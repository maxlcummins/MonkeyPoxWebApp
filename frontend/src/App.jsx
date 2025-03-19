// src/App.jsx
import { useCallback, useState } from "react";
import { useDropzone } from "react-dropzone";
import { uploadFiles } from "./api";

function App() {
  const [files, setFiles] = useState([]);
  const [message, setMessage] = useState("");
  const [uploading, setUploading] = useState(false);
  const [uploadProgress, setUploadProgress] = useState(0);

  const onDrop = useCallback((acceptedFiles, rejectedFiles) => {
    if (rejectedFiles.length > 0) {
      setMessage("Invalid file type. Only .fastq.gz files are allowed.");
      return;
    }
    setFiles(acceptedFiles);
    setMessage("");
    setUploadProgress(0);
  }, []);

  const handleUpload = async () => {
    if (files.length === 0) {
      setMessage("Please select files first.");
      return;
    }

    setUploading(true);
    try {
      const response = await uploadFiles(files, (progressEvent) => {
        const percentCompleted = Math.round(
          (progressEvent.loaded * 100) / progressEvent.total
        );
        setUploadProgress(percentCompleted);
      });

      console.log("Upload response:", response);
      setMessage(`Pipeline started! Run ID: ${response.run_id}`);
      setFiles([]);
    } catch (error) {
      console.error("Upload error:", error);
      if (error.message === "Network Error") {
        setMessage("Network Error");
      } else {
        setMessage(error.message || "Upload failed. Check your backend.");
      }
    } finally {
      setUploading(false);
      setUploadProgress(0);
    }
  };

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: ".fastq.gz", // Only allow .fastq.gz files
    maxSize: 1000000000, // 1GB limit (adjust as needed)
  });

  return (
    <div className="max-w-md mx-auto mt-10 p-4">
      <h1 className="text-3xl font-bold mb-6 text-center">FASTQ Upload</h1>

      <div
        {...getRootProps()}
        className={`border-2 border-dashed rounded p-6 text-center cursor-pointer ${
          isDragActive ? "bg-gray-200" : "hover:bg-gray-100"
        }`}
      >
        <input {...getInputProps()} />
        {isDragActive ? (
          <p>Drop the files here ...</p>
        ) : (
          <p>Drag & drop .fastq.gz files here, or click to select files</p>
        )}
      </div>

      {files.length > 0 && (
        <div className="mt-4">
          <h2 className="text-lg font-medium">Selected Files:</h2>
          <ul className="list-disc list-inside">
            {files.map((file, index) => (
              <li key={index}>{file.name}</li>
            ))}
          </ul>
        </div>
      )}

      {uploading && (
        <div className="mt-4">
          <p>Uploading... {uploadProgress}%</p>
          <progress value={uploadProgress} max="100" />
        </div>
      )}

      <button
        className="mt-6 w-full bg-blue-600 text-white py-2 rounded hover:bg-blue-700 disabled:opacity-50"
        onClick={handleUpload}
        disabled={uploading}
      >
        {uploading ? "Uploading..." : "Upload & Start Pipeline"}
      </button>

      {message && (
        <div className="mt-4 p-2 border rounded text-center">
          {message}
        </div>
      )}
    </div>
  );
}

export default App;