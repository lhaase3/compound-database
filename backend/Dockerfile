FROM python:3.10-slim

# Install system dependencies for RDKit drawing
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . /app

# Force NumPy <2.0 for RDKit compatibility before installing other requirements
RUN pip install "numpy<2.0"
RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8080

CMD ["gunicorn", "--bind", "0.0.0.0:8080", "app:app"]