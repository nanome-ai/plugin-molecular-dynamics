docker build -f test.Dockerfile -t test:latest ..

if [ "$(docker ps -aq -f name=test)" != "" ]; then
        # cleanup
        echo "removing exited container"
        docker rm -f test
fi

if [ "$1" != "" ]; then
    echo "Using specified plugin server: $1"
    docker run -d \
    -p 8888:8888 \
    -e PLUGIN_SERVER=$1 \
    --name test test
else
    echo "Using default plugin server: plugins.nanome.ai"
    docker run -d \
    -p 8888:8888 \
    --name test test
fi