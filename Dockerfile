FROM rust:1.82 as builder

WORKDIR /app

# Copy the Cargo.toml and Cargo.lock files to the container
COPY Cargo.toml Cargo.lock ./

# Copy the workspace members' Cargo.toml files
COPY definitions/Cargo.toml definitions/
RUN mkdir definitions/src && echo "fn test() {}" > definitions/src/lib.rs
COPY haplotyper/Cargo.toml ./haplotyper/
RUN mkdir haplotyper/src && echo "fn test_lib() {}" > haplotyper/src/lib.rs

COPY cli/Cargo.toml cli/
RUN mkdir cli/src && echo "fn main() {}" > cli/src/main.rs
COPY sandbox/Cargo.toml sandbox/
RUN mkdir sandbox/src && echo "fn main() {}" > sandbox/src/main.rs

# Build the dependencies to leverage Docker cache
RUN cargo build --release

# Copy the source code to the container
RUN rm -rf definitions/ haplotyper/ sandbox/ cli/
COPY . .
RUN touch definitions/src/lib.rs haplotyper/src/lib.rs 

RUN cargo build --release
# Stage 2: Create the runtime image
FROM debian:bookworm
COPY --from=builder /app/target/release/jtk /usr/bin/
RUN apt-get update && apt-get install -y curl unzip &&\
    curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" &&\
    unzip awscliv2.zip && ./aws/install && rm -rf ./aws
CMD ["jtk"]
