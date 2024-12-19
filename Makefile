.PHONY: \
	download-external-sources \
	docker-image \
	docker-interactive \
	docker-container \
	docker singularity

download-external-sources:
	./scripts/download-external-sources.sh

docker-image:
	@docker image build -t metaline:latest .

docker-container:
	@docker container run -dit --name metaline metaline:latest

docker-interactive:
	@docker image build -t metaline:latest .
	@docker container -dit --name metaline metaline:latest

singularity:
	./scripts/get_singularity_def_file_from_dockerfile.sh
	sudo singularity build metaline.sif metaline-singularity.def

