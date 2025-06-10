import type { NextConfig } from 'next'

const nextConfig: NextConfig = {
  async rewrites() {
    return [
      {
        source: '/indigo-api/:path*',
        destination: 'http://localhost:8002/:path*', // This proxies to the Docker Indigo service
      },
    ]
  },
  // Include any other config settings you're using
}

export default nextConfig
