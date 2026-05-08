import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    proxy: {
      // Forward all /api, /databases, /run, /results, /status calls to FastAPI backend
      '/databases': 'http://localhost:8000',
      '/run':       'http://localhost:8000',
      '/results':   'http://localhost:8000',
      '/status':    'http://localhost:8000',
      '/update':    'http://localhost:8000',
      '/license':   'http://localhost:8000',
    },
  },
  build: {
    outDir: 'dist',
  },
})
