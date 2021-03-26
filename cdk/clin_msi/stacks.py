from aws_cdk import core
from .codepipeline import ClinMsiCodePipeline


class CICDStack(core.NestedStack):

    def __init__(self, scope: core.Construct, id: str, **kwargs) -> None:
        super().__init__(scope, id, **kwargs)
        self.cicd = ClinMsiCodePipeline(self, 'Pipeline')


class ClinMsiStack(core.Stack):

    def __init__(self, scope: core.Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)
        
        # codepipeline
        self.cicd = ClinMsiCodePipeline(self, 'CICD', artifact_bucket_name="codepipeline-us-east-2-232898350506")

        # cf template output
        core.CfnOutput(self, 'OutputArtifactBucket', value=self.cicd.artifact_bucket.bucket_arn, description='Artifact Bucket ARN')
        core.CfnOutput(self, 'OutputPipelineArn', value=self.cicd.code_pipeline.pipeline_arn, description="CodePipeline ARN")
        core.CfnOutput(self, 'OutputPytestCodeBuildArn', value=self.cicd.pytest_project.project_arn, description="CodeBuild Pytest ARN")
        core.CfnOutput(self, 'OutputPublishCodeBuildArn', value=self.cicd.publish_project.project_arn, description="CodeBuild Publish ARN")
