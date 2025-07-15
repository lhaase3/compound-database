
# 1. Build stage
FROM node:18 AS builder

WORKDIR /app
COPY . .

# Install dependencies and build
RUN npm install --legacy-peer-deps
RUN npm run build

# Verify build output
RUN ls -l /app/.next

# 2. Run stage
FROM node:18-slim

# Set production environment
ENV NODE_ENV=production
ENV PORT=8080

WORKDIR /app

COPY --from=builder /app/package.json ./
COPY --from=builder /app/.next ./.next
COPY --from=builder /app/public ./public
COPY --from=builder /app/node_modules ./node_modules

EXPOSE 8080

CMD ["npm", "start"]
