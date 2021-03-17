from aws_cdk import core
from .codepipeline import ClinMsiCodePipeline


class CICDStack(core.NestedStack):

    def __init__(self, scope: core.Construct, id: str, **kwargs) -> None:
        super().__init__(scope, id, **kwargs)
        self.cicd = ClinMsiCodePipeline(self, 'Pipeline')


class ClinMsiStack(core.Stack):

    def __init__(self, scope: core.Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)
        stack = CICDStack(self, 'CICD')
        core.CfnOutput(self, 'OutputPipelineArn', value=stack.cicd.code_pipeline.pipeline_arn, description="CodePipeline ARN")
        core.CfnOutput(self, 'OutputPytestCodeBuildArn', value=stack.cicd.pytest_project.project_arn, description="CodeBuild Pytest ARN")
        core.CfnOutput(self, 'OutputPublishCodeBuildArn', value=stack.cicd.publish_project.project_arn, description="CodeBuild Publish ARN")
