.PHONY: \
	download-external-sources \
	docker-image \
	docker-interactive \
	docker-container \
	docker singularity

download-external-sources:
	./_scripts_singularity/download-external-sources.sh

docker-image:
	@docker image build -t metaline:latest .

docker-container:
	@docker container run -dit --name metaline metaline:latest

docker-interactive:
	@docker image build -t metaline:latest .
	@docker container -dit --name metaline metaline:latest

singularity:
	./_scripts_singularity/get_singularity_def_file_from_dockerfile.sh
	sudo singularity build metaline.sif metaline-singularity.def

