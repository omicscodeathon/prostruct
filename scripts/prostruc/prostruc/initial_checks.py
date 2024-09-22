import docker
from docker.errors import DockerException
import socket



def check_for_internet(host="8.8.8.8",port=53,timeout=3):
    try:
        socket.create_connection(address=(host,port),timeout=timeout)
        return True

    except OSError:
        return False


def check_for_docker():
    try:
        client = docker.from_env()
        version = client.version()
        # print(f"[*] Docker {version} is available")
        return True

    except DockerException as error:
        print(f"{error} Docker is either not available or not running")
        return False


def docker_list_containers():
    client = docker.from_env()
    containers = client.containers.list(all=True)
    if containers:
        for container in containers:
            print(f"Container ID: {container.id}, Image: {container.image.tags}, Status: {container.status}")

    else:
        print("No containers found")


def docker_check_image(image_tag):
    client = docker.from_env()
    images = client.images.list()
    found_image = False

    if images:
        for image in images:
            if image.tags and image_tag in image.tags:
                print(f"[*] {image_tag} docker image found")
                found_image = True
                return True

        if not found_image:
            print(f"ERROR: {image_tag} not found")
            print(f"[*] Pulling {image_tag}")
            docker_pull_image(image=image_tag)
    else:
        print("ERROR: No image found")
        print(f"[*] Pulling {image_tag}")
        docker_pull_image(image=image_tag)


def docker_pull_image(image):
    try:
        client = docker.from_env()
        client.images.pull(repository=image)
        print(f"Image: {image} pulled successfully")
        return True

    except Exception as error:
        print(f"ERROR: Failed to pull {image}")
        print(f"ERROR: {error}")
        return False


