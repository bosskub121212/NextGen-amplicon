# 16S/12S Amplicon App — Installation Guide
> Version: May 2025 | Tested on Ubuntu 22.04 (WSL2 + native Linux)

---

## สิ่งที่เครื่องปลายทางต้องมี

| รายการ | ขั้นต่ำ | แนะนำ |
|--------|---------|--------|
| RAM | 8 GB | 16 GB+ |
| Storage ว่าง | 20 GB | 50 GB+ |
| CPU | 4 cores | 8+ cores |
| OS | Ubuntu 20.04 / WSL2 | Ubuntu 22.04 / WSL2 |
| Internet | ต้องการ (ช่วง setup ครั้งแรก) | — |

> **Windows:** ต้องใช้ผ่าน **WSL2** เท่านั้น  
> **Linux / Mac:** รันได้โดยตรง

---

## ขั้นตอนที่ 0 — เตรียม WSL2 (Windows เท่านั้น)

เปิด **PowerShell** แบบ Administrator แล้วรัน:

```powershell
wsl --install
```

Restart เครื่อง จากนั้นเปิด **Ubuntu** จาก Start Menu แล้วตั้ง username/password

### ให้ WSL2 ใช้ RAM เพียงพอ (สำคัญมาก)

สร้างไฟล์ `C:\Users\<ชื่อผู้ใช้>\.wslconfig` แล้วใส่:

```ini
[wsl2]
memory=8GB
swap=16GB
processors=4
```

จากนั้นรันใน PowerShell:
```powershell
wsl --shutdown
```
แล้วเปิด Ubuntu ใหม่

---

## ขั้นตอนที่ 1 — แตก zip

เปิด Ubuntu (WSL2) แล้วรัน:

```bash
# แตก zip ไปที่ home directory
unzip /mnt/c/Users/<ชื่อ>/Downloads/amplicon_app_*.zip -d ~/amplicon_app

cd ~/amplicon_app
```

---

## ขั้นตอนที่ 2 — ติดตั้ง (ครั้งเดียว ใช้เวลา 20–40 นาที)

```bash
bash setup.sh
```

script จะติดตั้งอัตโนมัติ:
- Python 3 + FastAPI + uvicorn + psutil
- R 4.x + DADA2 (Bioconductor) + ggplot2 + vegan + ape + pheatmap
- Node.js 20 + build React frontend

---

## ขั้นตอนที่ 3 — ดาวน์โหลด Taxonomy Database

ต้องทำก่อนรันครั้งแรก — SILVA **ไม่ได้รวมอยู่ในไฟล์** เพราะขนาดใหญ่เกินไป

```bash
cd ~/amplicon_app/backend/databases

# สำหรับ 16S (แนะนำ) — ~130 MB
wget https://zenodo.org/records/10403693/files/silva_nr99_v138.2_train_set.fa.gz

# Species assignment (optional) — ~78 MB
wget https://zenodo.org/records/10403693/files/silva_species_assignment_v138.2.fa.gz
```

สำหรับ 12S / ITS / database อื่น: วางไฟล์ `.fa.gz` หรือ `.fasta.gz` ในโฟลเดอร์ `backend/databases/` แล้ว app จะ detect อัตโนมัติ

---

## ขั้นตอนที่ 4 — เริ่มใช้งาน

```bash
cd ~/amplicon_app
bash start.sh
```

เปิด browser แล้วไปที่: **http://localhost:5173**

หยุด app: กด `Ctrl+C`

---

## การใช้งานทุกครั้ง

```bash
cd ~/amplicon_app
bash start.sh
```

---

## โครงสร้างไฟล์

```
amplicon_app/
├── backend/
│   ├── main.py                  # FastAPI backend
│   ├── requirements.txt         # Python packages
│   ├── r_scripts/
│   │   ├── dada2_pipeline.R     # DADA2 pipeline หลัก
│   │   └── replot.R             # สร้างกราฟใหม่ด้วยสีที่กำหนด
│   └── databases/               # ← วาง SILVA / database ที่นี่
├── frontend/
│   ├── src/                     # React source code
│   └── dist/                    # Built frontend (พร้อมใช้งาน)
├── venv/                        # Python venv (สร้างตอน setup)
├── setup.sh                     # ติดตั้งทุกอย่าง
├── start.sh                     # เริ่มใช้งาน
├── start_dev.sh                 # สำหรับนักพัฒนา (hot-reload)
└── INSTALL.md                   # คู่มือนี้
```

---

## แก้ปัญหาที่พบบ่อย

### ❌ Memory error ตอน taxonomy assignment
SILVA ต้องการ RAM ~5–8 GB เพิ่ม swap:
```bash
sudo fallocate -l 16G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
# ทำให้ถาวร:
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

### ❌ R package install ล้มเหลว
```bash
# ติดตั้ง R packages ใหม่ทีละตัว
Rscript -e "options(repos='https://cloud.r-project.org'); install.packages('ggplot2')"
Rscript -e "BiocManager::install('dada2')"

# หรือรัน script ครบชุด
Rscript ~/amplicon_app/install_r_packages.R
```

### ❌ Port ชน / app ไม่ start
```bash
pkill -f uvicorn
pkill -f "serve dist"
sleep 2
bash ~/amplicon_app/start.sh
```

### ❌ Backend ขึ้น error อื่น
```bash
# ดู log โดยตรง
cd ~/amplicon_app/backend
source ../venv/bin/activate
uvicorn main:app --host 0.0.0.0 --port 8000
```

### ❌ setup.sh ค้างหรือ timeout
ลอง restart WSL แล้วรัน setup.sh ใหม่ — มันจะ skip packages ที่ติดตั้งแล้ว

---

## หมายเหตุสำหรับ Linux / Mac

ขั้นตอนเหมือนกันทุกอย่าง แต่ข้ามขั้นตอน WSL2 ได้เลย  
สำหรับ Mac ให้ใช้ `brew install r` แทนการใช้ `apt`
