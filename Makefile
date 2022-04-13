

up:
	export USER_ID=$$(id -u $$USER); export USER_GROUP=$$(id -g $$USER); docker-compose up -d --build

down:
	docker-compose down
