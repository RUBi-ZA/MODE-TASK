pipeline {
    agent { docker { image 'python:3.5.1' } }
    options { timestamps(); timeout(time: 20, unit: 'MINUTES') }
    stages {
        stage('build') {
            steps {
                sh 'python --version'
            }
        }
    }
}
