import os

from aws_cdk import core
from aws_cdk import aws_iam as iam
from aws_cdk import aws_codebuild as cb
from aws_cdk import aws_codepipeline as cp
from aws_cdk import aws_codepipeline_actions as actions


class ClinMsiCodePipeline(core.Construct):
    """Defines a CodePipeline: GitHub -> CodeBuild (PyTest) -> Manual Approval -> CodeBuild (Deploys to PyPI)"""

    @property
    def pytest_project(self):
        return self._pytest_project
    
    @property
    def publish_project(self):
        return self._publish_project

    @property
    def code_pipeline(self):
        return self._pipeline

    def __init__(self, scope: core.Construct, id: str, **kwargs) -> None:
        """Define the resources for the CodePipeline.
        
        :param scope: the parent construct
        :param id: the logical id
        :param artifact_bucket_name: the bucket name for artifacts passed between code pipeline 
        """
        super().__init__(scope, id, **kwargs)

        account_id = core.Stack.of(self).account
        env = core.Stack.of(self).environment

        repo_name = 'clin-msi'

        # build projects
        environment = cb.BuildEnvironment(build_image=cb.LinuxBuildImage.STANDARD_4_0)
        self._pytest_project = cb.PipelineProject(self, 'UnitTest', 
            project_name=f"{repo_name}-unittest",
            build_spec=cb.BuildSpec.from_object(buildspec_pytest),
            environment=environment)
        self._publish_project = cb.PipelineProject(self, 'Publish',
            project_name=f"{repo_name}-publish",
            build_spec=cb.BuildSpec.from_object(buildspec_publish),
            environment=environment)
        self._publish_project.add_to_role_policy(
            iam.PolicyStatement(
                effect=iam.Effect.ALLOW,
                actions=["ssm:GetParameters"],
                resources=[
                    f"arn:aws:ssm:{env}:{account_id}:parameter/ClinMsi/PyPI/Credentials/Username",
                    f"arn:aws:ssm:{env}:{account_id}:parameter/ClinMsi/PyPI/Credentials/Password"
                ]
            ))

        # actions
        source_output = cp.Artifact()
        source_action = actions.BitBucketSourceAction(
            action_name='GitHub',
            connection_arn=f"arn:aws:codestar-connections:{env}:{account_id}:connection/a859d7f0-0bf3-48f5-bce1-492c9bc08bef",
            output=source_output,
            code_build_clone_output=True,
            owner="nch-igm",
            repo=repo_name,
            branch="master")

        unittest_action = actions.CodeBuildAction(action_name='UnittestAction', input=source_output, project=self._pytest_project)
        publish_action = actions.CodeBuildAction(action_name='PyPIPublishAction', input=source_output, project=self._publish_project)
        manual_approval = actions.ManualApprovalAction(action_name='ApproveToPublish')

        # initialize pipeline
        self._pipeline = cp.Pipeline(self, 'CodePipeline', pipeline_name=repo_name)
        self._pipeline.add_to_role_policy(
            iam.PolicyStatement(
                effect=iam.Effect.ALLOW,
                actions=['codebuild:StartBuild', 'kms:PutKeyPolicy'],
                resources=[
                    self._pytest_project.project_arn,
                    self._publish_project.project_arn,
                    self._pipeline.artifact_bucket.bucket_arn
                ]))
        self._pipeline.add_stage(stage_name="Source", actions=[source_action])
        self._pipeline.add_stage(stage_name="UnitTest", actions=[unittest_action])
        self._pipeline.add_stage(stage_name="ManualApprovalForPublish", actions=[manual_approval])
        self._pipeline.add_stage(stage_name="Publish", actions=[publish_action])

    
buildspec_pytest = {
    "version": '0.2',
    "phases": {
        "install": {
            "runtime-versions": {
                "python": "3.8",
                "nodejs": "10"
            },
            "commands": [
                "npm install -g aws-cdk",
                "pip3 install poetry"
            ]
        },
        "build": {
            "commands": [
                "poetry install -vvv",
                "poetry run pytest tests/"
            ]
        },
    }
}

buildspec_publish = {
    "version": '0.2',
    "env": {
        "parameter-store": {
            "PYPI_USERNAME": "/ClinMsi/PyPI/Credentials/Username",
            "PYPI_PASSWORD": "/ClinMsi/PyPI/Credentials/Password" 
        }
    },
    "phases": {
        "install": {
            "runtime-versions": {
                "python": "3.8",
                "nodejs": "10"
            },
            "commands": [
                "npm install -g aws-cdk",
                "pip3 install poetry"
            ]
        },
        "build": {
            "commands": [
                "poetry build",
                "poetry publish -u $PYPI_USERNAME -p $PYPI_PASSWORD"
            ]
        },
    }
}